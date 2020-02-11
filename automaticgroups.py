import sympy
import numpy as np
import subprocess32
import subprocess
import os
import time
import glob
import itertools
import copy
import random

def certify_hyperbolicity(relator,tryhard=1,generators=None,timeout=20,verbose=False,cleanup=True,**kwargs):
    """
    Attempt to check if a one relator group is hyperbolic. 
    If return is True then group is hyperbolic. If return is False then test was inconclusive, try longer timeout.
    'relator' is a string whose characters all belong to 'generators'. 
    'generators' is list of letters, closed under case change. Ordering defines lexicographic order. If not specified then will be generated in numerical order, ie {'B','A,'a','b']
    If 'cleanup=True' then delete all the files created for/by kbmag.
    First attempt to find a shortlex automatic structure, then attempt to check hyperbolicity.
    In principle, if the group is hyperbolic both of these steps will succeed, given enough time, and the function will return True.
    However, if the group is shortlex automatic but not hyperbolic the second step will not terminate.
    'timeout' is time in seconds to allow each of the two steps to run.
    tryhard=0 just tries given generators and timeout and quits if unsuccessful. 
    tryhard=1 tries once, and if inconclusive tries again with double timeout time
    tryhard=2 tries once, and if inconclusive tries 10 more times with double timeout time and random permutation of generator order
    tryhard=3 tries once, and if inconclusive tries again with double timeout time and every possible permutaiton of generator order
    """
    if type(relator)==str:
        relatorasstring=relator
        relatoraslist=letterstringtointlist(relator)
    elif type(relator)==list and type(relator[0])==int:
        relatoraslist=relator
        relatorasstring=intlisttoletterstring(relator)
    else:
        raise UsageError('relator should be a string or list of nonzero integers')
    if not generators:
        rank=max(abs(x) for x in relatoraslist)
        generators=[intlisttoletterstring([i]) for i in range(-rank,0)+range(1,rank+1)]
    if 'tmp_directory' in kwargs:
        directory=kwargs['tmp_directory']
    else:
        directory='kbmag_tmp_directory'
    Icreatedthisdirectory=False
    if not os.path.exists(directory):
        os.mkdir(directory)
        Icreatedthisdirectory=True
    if 'filename' in kwargs:
        thefilename=kwargs['filename']
    else:
        thefilename="OneRelatorGroup-"+relatorasstring+"-"+"".join(generators)
    writetokbmagfile(directory+'/'+thefilename,generators,[relatorasstring])
    aut=False 
    hyp=False 
    if verbose:
        print "Attempting to find automatic structure with generator order: "+str(generators)
    try:
        #subprocess32.run(['autgroup','-silent',directory+'/'+thefilename],check=True,timeout=timeout)
        if 'kbprogargs' in kwargs:
            kbprogargument=['kbprog']+kwargs['kbprogargs']+[directory+'/'+thefilename]
        else:
            kbprogargument=['kbprog','-mt', '20', '-hf', '100', '-cn','0', '-me', '200000', '-silent', '-wd',directory+'/'+thefilename]
        subprocess32.run(kbprogargument,check=True,timeout=timeout)
        subprocess32.run(['gpmakefsa','-silent', directory+'/'+thefilename],check=True,timeout=timeout)
        subprocess32.run(['gpaxioms','-silent',directory+'/'+thefilename],check=True,timeout=timeout)
        aut=True # all subprocesses completed in time and with returncode=0
        if verbose:
            print "Automatic structure found."
    except (subprocess32.TimeoutExpired,subprocess32.CalledProcessError) as e: # if either autgroup timed out or complete with nonzero returncode
        if verbose:
            print "Failed to find automatic structure with error: "+str(e)
        pass # aut remains False
    if aut:
        try:
            if verbose:
                print "Checking hyperbolicity."
            subprocess32.run(['gpgeowa','-silent',directory+'/'+thefilename],check=True,timeout=timeout)
            hyp=True
        except (subprocess32.TimeoutExpired,subprocess32.CalledProcessError) as e:
            if verbose:
                print "Failed to find hyperbolic structure with error: "+str(e)
            pass # hyp remains false
    if cleanup:
        files = glob.glob(directory+'/'+thefilename+"*")
        for file in files:
            try:
                os.remove(file)
            except:
                pass
        if Icreatedthisdirectory:
            try:
                os.rmdir(directory)
            except OSError:# sometimes files escape the cleanup. I think this happens because subprocess has timed out and python moves on and tries to cleanup before kbprog finishes and saves its files. Wait and try again.
                time.sleep(5)
                files = glob.glob(directory+'/'+thefilename+"*")
                for file in files:
                    try:
                        os.remove(file)
                    except:
                        pass
                if Icreatedthisdirectory:
                    try:
                        os.rmdir(directory)
                    except:
                        pass
    if aut and hyp:
        return True
    else:
        if tryhard==1:
            if verbose:
                print "Trying again with double wait time."
            return certify_hyperbolicity(relator,0,generators,2*timeout,verbose,cleanup,**kwargs)
        if tryhard==2:
            if verbose==True:
                for i in range(10):
                    orderedgens=random.sample(generators,len(generators))
                    print "Trying kbmag with generator order "+str(orderedgens)
                    result=certify_hyperbolicity(relator,0,orderedgens,2*timeout,verbose,cleanup,**kwargs)
                    if result==True:
                        return True
                else:
                    return False
            else:
                return any(certify_hyperbolicity(relator,0,random.sample(generators,len(generators)),2*timeout,verbose,cleanup) for i in range(10))
        if tryhard==3:
            return any(certify_hyperbolicity(relator,0,permutedgens,2*timeout,verbose,cleanup) for permutedgens in itertools.permutations(generators))
        return False



def growthseries(generators,relator,verbose=False,cleanup=True,**kwargs):
    """
    Compute the growth series of the shortlex automatic one relator group with given generators and relator.
    'generators' is list of letters, closed under case change. Ordering defines lexicographic order.
    'relator' is a string whose characters all belong to 'generators'. 
    If 'cleanup=True' then delete all the files created for/by kbmag.
    Output is two lists of integers that are coefficients of numerator and denominator of rational function whose Taylor expansion is the growth series. Lists are ordered with higher order coefficients first.
    """
    if 'filename' in kwargs:
        thefilename=kwargs['filename']
    else:
        thefilename="OneRelatorGroup-"+relator
    writetokbmagfile(thefilename,generators,[relator])
    if verbose:
        autrun=subprocess.call(['autgroup','-v',thefilename])
    else:
        autrun=subprocess.call(['autgroup','-silent',thefilename])
    assert(autrun==0)
    if verbose:
        grrun=subprocess.call(['fsagrowth','-v',thefilename+'.wa'])
    else:
        grrun=subprocess.call(['fsagrowth',thefilename+'.wa'])
    assert(grrun==0)
    rationalfunction=fsagrowthtopolystrings(thefilename+'.wa.growth')
    if cleanup:
        files = glob.glob('./'+thefilename+"*")
        for file in files:
            os.remove(file)
    return rationalfunction

def growthrate(generators,relator,verbose=False,cleanup=True,**kwargs):
    """
    Compute the growth exponent of the shortlex automatic one relator group with given generators and relator.
    'generators' is list of letters, closed under case change. Ordering defines lexicographic order.
    'relator' is a string whose characters all belong to 'generators'. 
    If 'cleanup=True' then delete all the files created for/by kbmag.
    This function reads the rational function encoding the growth series and returns the reciprocal of the smallest positive real root of the denominator. 
    """
    return 1/smallpole(*growthseries(generators,relator,verbose,cleanup,**kwargs))

def numericalgrowthrate(generators=None,relator=None,verbose=False,cleanup=True,**kwargs):
    """
    Compute the growth exponent of the shortlex automatic one relator group with given generators and relator.
    'generators' is list of letters, closed under case change. Ordering defines lexicographic order.
    'relator' is a string whose characters all belong to 'generators'. 
    If 'cleanup=True' then delete all the files created for/by kbmag.
    This function gets the (transpose of) the transition table of the shortlex automaton and calculates the largest real eigenvalue numerically.
    """
    return largestrealeigenvalue(automatatransitionmatrix(generators,relator,verbose,cleanup,**kwargs))


def automatatransitionmatrix(generators=None,relator=None,verbose=False,cleanup=True,**kwargs):
    """
    Return (transpose of the) transition matrix of the shortlex automata for one relator group with given generators and relator.
    'generators' is list of letters, closed under case change. Ordering defines lexicographic order.
    'relator' is a string whose characters all belong to 'generators'. 
    If 'cleanup=True' then delete all the files created for/by kbmag.
    if 'inputwafile' is specified then generators and relators are ignored and the transition table is read directly from an existing file. Otherwise such a .wa file is created by running autgroup from kbmag.
    """
    if 'inputwafile' in kwargs:
        f=open(kwargs['inputwafile'],'r')
    else:
        if 'filename' in kwargs:
            thefilename=kwargs['filename']
        else:
            thefilename="OneRelatorGroup-"+relator
        writetokbmagfile(thefilename,generators,[relator])
        if verbose:
            autrun=subprocess.call(['autgroup','-v',thefilename])
        else:
            autrun=subprocess.call(['autgroup','-silent',thefilename])
        assert(autrun==0)
        f=open(thefilename+".wa","r")
    lines=f.readlines()
    currentline=0
    for thisline in lines:
        if "states :=" in thisline:
            break
        else:
            currentline=currentline+1
    for i in range(currentline,len(lines)):
        if "size :=" in lines[i]:
            size=int(lines[i].split("size := ",1)[1].rstrip(' ,\n'))
            break
    else:
        raise RunTimeError("Couln't find size of the automata.")
    currentline=0
    for thisline in lines:
        if " table := rec(" in thisline:
            break
        else:
            currentline=currentline+1
    else:
        raise RunTimeError("Couldn't find the table record.")
    for i in range(currentline,len(lines)):
        thisline=lines[i]
        if "format := " in thisline:
            if "dense deterministic" in thisline:
                theformat='dense'
            elif "sparse" in thisline:
                theformat='sparse'
                NotImplemented
            else:
                raise RunTimeError("Couldn't find format string.")
            break
    else:
        raise RunTimeError("Couldn't find format string.")
    currentline=0
    for thisline in lines:
        if "transitions := [" in thisline:
            break
        else:
            currentline=currentline+1
    else:
        raise RunTimeError("Couldn't find the transition table.")
    transitiontable=''
    openparens=1
    transtableline=lines[currentline].split('transitions := [',1)[1]
    while openparens:
        for x in transtableline:
            if x=='[':
                openparens=openparens+1
            elif x==']':
                openparens=openparens-1
                if openparens==0:
                    break
            if x==' ' or x=='\n' or x=='\t':
                pass
            elif x==',' and transitiontable[-1]==']':
                transitiontable+=';'
            else:
                transitiontable+=x
        if openparens:
            currentline=currentline+1
            transtableline=lines[currentline]
    trows=transitiontable.split(';')
    assert(size==len(trows))
    tmat=[]
    for trow in trows:
        targets=[int(c)-1 for c in  (trow.strip(' ,[]\n')).split(',') if c!='' and int(c)>0]
        tmatrow=[]
        for j in range(size):
            if j in targets:
                tmatrow.append(1)
            else:
                tmatrow.append(0)
        tmat.append(tmatrow)
    if cleanup and not 'inputwafile' in kwargs:
        files = glob.glob('./'+thefilename+"*")
        for file in files:
            os.remove(file)
    return tmat








def wordreduce(thestring,thefilename):
    """
    Given thestring returns the shortlex minimal string representing the same group element as determined by automatic strucure for the group defined in thefilename.
    """
    p=subprocess.Popen(['wordreduce',thefilename],stdin=subprocess.PIPE,stderr = subprocess.STDOUT,stdout=subprocess.PIPE)
    output=p.communicate(addstars(thestring)+';')
    return removestars((output[0].rstrip('\n'))[95:])


class groupelement(object):
    """
    Defines an element in an automatric group defined in the kbmag file 'thefilename'.
    Element represented by self.string which is shortlex minimal representation of that group element as determined by the automatic structure.
    """
    def __init__(self,thestring,thefilename):
        if os.path.isfile(thefilename):
            self.groupfilename=thefilename
        else:
            raise NameError('group file not found')
        if not os.path.isfile(thefilename+'.diff1'):
            ec=subprocess.call(['autgroup','-silent',thefilename])
            if ec!=0:
                raise RuntimeError('autgroup failed with value '+str(ec))
        try:
            self.string=wordreduce(thestring,thefilename)
        except OSError:#sometimes wordreduce fails for no apparent reason and it is sufficient to just try again
            try:
                self.string=wordreduce(thestring,thefilename)
            except OSError:
                try:
                    self.string=wordreduce(thestring,thefilename)
                except OSError:
                    print thestring,thefilename
                    assert(False)
                
    def __str__(self):
        return self.string

    def __repr__(self):
        return self.string

    def inverse(self):
        return groupelement(''.join([x for x in reversed(self.string.swapcase())]),self.groupfilename)

    def __pow__(self,thepower):
        if thepower==0:
            return groupelement('',self.groupfilename)
        elif thepower<0:
            return (self.inverse())**(-thepower)
        else:
            return groupelement(self.string*thepower,self.groupfilename)
        
    def __mul__(self,other):
        return groupelement(self.string+other.string,self.groupfilename)
        
    def __eq__(self,other):
        return self.string==other.string


#------------------- Auxiliary functions for interacting with kbmag
    
def writetokbmagfile(filename,generators,relators,**kwargs):
    """
    Create a kbmag file containing group definition.
    'generators' should be a list of generators consisting of lowercase letters and their inverses the corresponding uppercase letters. The order in which they are given is used by kbmag to define lexicographic ordering on words. 
    'relators' should be list of strings conisisitng of upper and lower case letters corresponding to generators and their inverses.
    """
    f=open(filename,"w")
    if 'description' in kwargs:
        f.write('# '+kwargs['description']+'\n')
    else:
        f.write('# some group\n')
    f.write('_RWS := rec(\n')
    f.write('  isRWS := true,\n')
    f.write('  ordering := "shortlex",\n')
    f.write('  generatorOrder := ['+','.join(generators)+'],\n')
    f.write('  inverses := ['+','.join([x.swapcase() for x in generators])+'],\n')
    f.write('  equations := ['+','.join(['['+str(addstars(relatorstring))+',IdWord]' for relatorstring in relators])+']\n')
    f.write(');')

def addstars(s):
    if s=='' or s==[]:
        return 'IdWord'
    else:
        t=''
        for i in range(len(s)-1):
            t=t+s[i]+'*'
        t=t+s[-1]
        return t

def removestars(s):
    if s=='IdWord':
        return ''
    else:
        rstring=s
        while '^' in rstring:
            i=rstring.index('^')
            theletter=rstring[i-1]
            if '*' not in rstring[i+1:]:
                thepower=int(rstring[i+1:])
                j=len(rstring)
            else:
                j=i+1+rstring[i+1:].index('*')
                thepower=int(rstring[i+1:j])
            rstring=rstring[:i-1]+theletter*thepower+rstring[j+1:]
        string=''
        for x in [x for x in rstring if x!='*']:
            string=string+x
        return string


def intlisttoletterstring(intlist):
    alphabet=['Z','Y','X','W','V','U','T','S','R','Q','P','O','N','M','L','K','J','I','H','G','F','E','D','C','B','A','','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    return ''.join([alphabet[26+x]  for x in intlist])

def letterstringtointlist(thestring):
    alphabet=['Z','Y','X','W','V','U','T','S','R','Q','P','O','N','M','L','K','J','I','H','G','F','E','D','C','B','A','','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    return [alphabet.index(c)-26 for c in thestring]

def fsagrowthtopolystrings(filename):
    f=open(filename,"r")
    lines=f.readlines()
    assert('as previous' in lines[-1] and 'as previous' in lines[-2])
    liness=lines[3:-2]
    line=''
    for l in liness:
        line=line+(l.rstrip('\n'))
    numbegin=line.index('(')
    numend=line.index(')')
    denombegin=1+numend+(line[1+numend:]).index('(')
    denomend=1+numend+(line[1+numend:]).index(')')
    return sympy.poly(line[1+numbegin:numend]).all_coeffs(),sympy.poly(line[1+denombegin:denomend]).all_coeffs()

def smallpole(num,denom,force=False):
    droots=np.roots(denom)
    smallrealdenomroot=droots[(droots.imag==0)&(droots.real>0)].real.min()
    if (not force) and np.isclose(np.polyval(np.array(num).astype(float),smallrealdenomroot),0):
        raise InputError('It appears that the smallest root '+str(smallrealdenomroot)+' of the denominator is also a root of the numerator.')
    return smallrealdenomroot




def largestrealeigenvalue(M):
    eigvals=np.linalg.eigvals(M)
    realeigvals=[x.real for x in eigvals if np.isreal(x)]
    return max(realeigvals)
