#!/usr/bin/env python
"""
Created on Mon 9 Sep 2013

author: Kristoffer Braekken
"""
import os,sys,time

"""CONSTANTS"""

_HEADER_EXTENSIONS = ['.h','.hpp']

"""CLASSES"""

class Project(object):
    """
    Contains the metainformation from the project.
    """
    def __init__(self, metafile):
        """
        @param metafile The file containing metainformation for the project.
        """
        self.data = {}
        self._readData(metafile)

    def __call__(self):
        """
        Can get hold of dictionairy by calling. This results in shorter code
        usage.
        """
        return self.data

    def _readData(self, metafile):
        """
        Stores the data in the instance. Private method.

        @param metafile File containing the data.
        """
        reader = open(metafile,'r')
        reader.readline() # Skip header
        for line in reader:
            if len(line.strip().split(':')) == 2:
                if line.strip().split(':')[0] in ['files','compileFlags','linkFlags',
                        'libLocations','libs','headerSearchDirs']:
                    # Some params are lists
                    identity = line.strip().split(':')[0]
                    value = line.strip().split(':')[1]
                    value = [i.strip() for i in value.strip('[').strip(']').split(',')]
                    # If only contains empty string, make it empty
                    if value[0] == '' and value[-1] == '':
                        value = []
                else:
                    identity,value = line.strip().split(':')
                self.data[identity] = value
            else:
                print 'Found line longer or shorter than one separation.'
                print 'Please correct metafile.'
                reader.close()
                sys.exit(1)

        reader.close()

if __name__ == '__main__':
    # Clear screen first
    os.system('clear')
    
    # Check that metafile is present
    for contents in os.listdir('.'):
        if contents == 'meta.dat':
            break;
        if contents == os.listdir('.')[-1] and contents != 'meta.dat':
            print 'No meta.dat file found in root of project (or here).'
            sys.exit(1)
    # Load metafile
    p = Project('meta.dat')

    # Fetch compiler
    compiler = p()['compiler']
    
    # Check if build dir present, if not, create
    if not os.path.exists(p()['buildDir']):
        os.mkdir(p()['buildDir'])
    
    # Check if objects dir present, if not create
    # --> At this stage, build dir will be present
    if not os.path.exists(os.path.join(p()['buildDir'],p()['objectsDir'])):
        os.mkdir(os.path.join(p()['buildDir'],p()['objectsDir']))
    
    start = time.time()
    
    # Create link and compile flagstrings
    cmpFlagString = ''; lnkFlagString = ''
    for flag in p()['compileFlags']:
        cmpFlagString += '-%s' % flag
        if not flag == p()['compileFlags'][-1]:
            cmpFlagString += ' '
    for path in p()['headerSearchDirs']:
        cmpFlagString += '-I%s' % path
        if not path == p()['headerSearchDirs'][-1]:
            cmpFlagString += ' '
    for flag in p()['linkFlags']:
        lnkFlagString += '-%s' % flag
        if not flag == p()['linkFlags'][-1]:
            lnkFlagString += ' '
    
    # Check if executable already present
    for existing_file in os.listdir(p()['buildDir']):
        if existing_file == '%s.x' % p()['exe']:
            # If it is, delete it
            os.remove('%s%s.x' % (p()['buildDir'],p()['exe']))
            print 'Deleted old executable.\n'
    
    print 'Compiling files:\n______________________','\n'
    
    # Compile
    for cfile in p()['files']:
        name = cfile.split('.')[0]
    
        if cmpFlagString != '':
            runstring = '%s %s -c %s%s -o %s%s.o' % (compiler, cmpFlagString,p()['sourceDir'],
                                                    cfile,os.path.join(p()['buildDir'],p()['objectsDir']),name)
        else:
            runstring = '%s -c %s%s -o %s%s.o' % \
                                (compiler, p()['sourceDir'],cfile,os.path.join(p()['buildDir'],p()['objectsDir']),name)
    
        print '$', runstring
    
        """
        Only compile this file if:
        --> There doesn't exist an object file with that name
        --> If the object file that exists has a creation date older than the
        modification date of the file to be compiled.
        """
        doCompile = False
    
        # File exists?
        if (name + '.o') in os.listdir(os.path.join(p()['buildDir'],p()['objectsDir'])):
            # Creation time older than modification time of to be compiled?
            if os.path.getmtime(os.path.join(p()['sourceDir'],cfile)) > \
                                os.path.getmtime(os.path.join(p()['buildDir'],p()['objectsDir']+name+'.o')):
                doCompile = True
            """
            Header file might have changed even though the file itself hasnt.

            Only check if condition isnt satisfied above.
            """
            if not doCompile:
                for ex in _HEADER_EXTENSIONS:
                    if (name + ex) in os.listdir(p()['sourceDir']):
                        # If has header file, same check as with main file
                        if os.path.getmtime(os.path.join(p()['sourceDir'],name+ex)) > \
                                os.path.getmtime(os.path.join(p()['buildDir'],p()['objectsDir']+name+'.o')):
                            doCompile = True
            
            if doCompile:
                # If exists yet still needs to compile, delete old one
                os.remove(os.path.join(p()['buildDir'],p()['objectsDir'] + name + '.o'))
            else:
                print '--> No changes, continuing.'
        else:
            # If no object file present, compile
            doCompile = True
    
        if doCompile: os.system(runstring)
        print '\n'
    
    """LINKING"""
    
    # Flags
    if len(p()['linkFlags']) != 0:
        linkstring = '%s %s -o %s%s ' % (compiler, lnkFlagString,p()['buildDir'],p()['exe'])
    else:
        linkstring = '%s -o %s%s ' % (compiler, p()['buildDir'],p()['exe'])
    
    linkstring += '%s*.o ' % os.path.join(p()['buildDir'],p()['objectsDir'])
    
    # Extra lib locations
    if len(p()['libLocations']) != 0:
        for libPath in p()['libLocations']:
            linkstring += '-L%s ' % libPath
    
    if len(p()['libs']) != 0:
        for lib in p()['libs']:
            if lib == p()['libs'][-1]:
                linkstring += '-l%s' % lib
                continue
            linkstring += '-l%s ' % lib
    
    print '______________________\n\nNow linking:\n','$',linkstring,'\n'
    os.system(linkstring)
    print ''
    
    end = time.time()
    
    print '______________________\n\n--> Sucessfully compiled in %g s.\n______________________' % (end - start)
