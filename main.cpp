/*
   progressive BKZ library by NICT security fundemental lab.
 */


//comment out when publishing
//#define __develop
//#define __ide
#ifdef __develop
    #include <develroutines/precodes.cpp>
#endif

#include <lattice/pbkz.hpp>
#include <lattice/genlattice.cpp>


#ifdef __develop
    #include <extendlib/extend.hpp>
    #include <develroutines/devel.hpp>
#endif
#include "procs.cpp"


int main(int argc, char** argv) {

    int i;
    std::vector<std::string> args;
    
#ifndef __ide
    //Extracting command-line options
    if (argc<=1) {
        cout << "usage: ./a.out -if [input file name] -of [output file name]" << endl;
        cout << "       -com [command] " << endl;
        return -1;
    }
    
    //Copying commands
    for (i=0;i<argc;i++) {
        std::string a;
        a = argv[i];
        args.push_back(a);
    }
#endif

#ifdef __ide
    std::string cl;
    //Simulate
    
    std::vector<std::string> args2;
    std::string t;
    char qflag=0;
    boost::algorithm::split(args2,cl,boost::is_any_of(" "));
    //Process double-quotations
    for (i=0;i<args2.size();i++) {
        if ((args2[i][0]=='\"') && (qflag==0)){
            t = args2[i].substr(1,args2[i].length());
            qflag = 1;
        } else {
            if (qflag==0) {
                args.push_back(args2[i]);
            } else {
                t += " " + args2[i];
                if (args2[i][args2[i].length()-1]=='\"') {
                    //last character is the double-quortation
                    args.push_back(t.substr(0,t.length()-1));
                    t = "";
                    qflag = 0;
                }
            }
        }
    }
    
    for (i=0;i<args.size();i++) {
        while (args[i].length()>0) {
            if (args[i][0]==' ') {
                args[i] = args[i].substr(1,args[i].length()-1);
            } else 
            if (args[i][args[i].length()-1]==' ') {
                args[i] = args[i].substr(0,args[i].length()-1);
            } else {
                break;
            }
        }
    }
    
#endif    
    
    for (i=0;i<(int)args.size();i++) {
        if ((args[i]=="-genbasis") || (args[i]=="-genmatrix") || (args[i]=="-lll") || (args[i]=="-bkz") || (args[i]=="-genstrategy") || (args[i]=="-enum") || (args[i]=="-pfunc") || (args[i]=="-disp") || (args[i]=="-randomizebasis") || (args[i]=="-simgslength")) {
            args[i] = "com=" + args[i].substr(1,args[i].length());
        }
        #ifdef __develop
            args[i] = othercommands(args[i]); 
        #endif
    }
    processcommand(args);

}
 
