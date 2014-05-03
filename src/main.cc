#include "kongrad.hh"
#define DEFINE_GLOBAL
#include "global.h"
#include <iostream>
#include <cassert>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

/** 
 * 
 * 
 * @mainpage CP 3: Konjugierte Gradienten
 * Lösung von linearen Gleichungssystemen mit der Methode der konjugierten Gradienten. Dieses Programm soll entsprechend der CP 3 Vorlesung 2014 geschrieben werden.
 *   @par Über die Vorlesung hinausgehende Ansprüche:
 *   - unit tests für möglichst alle Funktionen/Funktionalitäten
 *   - leserlicher Code
 *   - gute Dokumentation
 *   - Versionsverwaltung
 *
 * The repository can be found at https://github.com/micgro42/KonGrad
 * 
 *   @author Michael Große
 *
 * The website for the course can be found at: http://www-com.physik.hu-berlin.de/~bunk/cp3/
 * 
 */

/**
 * 
 * @file main.cc
 * 
 * @brief contains the main routine to call/test the conjugate gradient algorithm.
 * 
 * 
 * 
 */


namespace logging = boost::log;
using namespace std;

void init(const int loglevel){
    /// @todo so umschreiben, dass der Aufruf mit ./main --loglevel=info erfolgt
        switch (loglevel){
            case 0:
                logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::fatal);
                break;
            case 1:
                logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::error);
                break;
            case 2:
                logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::warning);
                break;
            case 3:
                logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::info);
                break;
            case 4:
                logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::debug);
                break;
            case 5:
                logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::trace);
                break;
            default:
                BOOST_LOG_TRIVIAL(fatal) << "first parameter must be int 0..5 or empty";
                assert(false);
                break;
        }
}

int main(int argc, char** argv){
    if (argc > 1){
        init(atoi(argv[1]));
    }else{
        init(3);
    }
    

    vector< vector<double> > A;
    vector<double> line;
    for (int i=0;i<3;++i){
        line.assign(3,0);
        line.at(i)=i+1;
        A.push_back(line);
    }
    
    vector<double> b;
    b.push_back(1);
    b.push_back(2);
    b.push_back(3);
    
    KonGrad LGS01(A, b);
    LGS01.testmv(b);
    LGS01.startRandomGenerator(0);
    
    
    
    vector<double> c;
    c.push_back(2);
    c.push_back(2);
    c.push_back(2);
    vector<double> result;
    LGS01.solve(c, result);
    


    ndim=2;
    //int * lsize = new int[3];
    lsize[1]=3;
    lsize[2]=3;




    return 0;
}


