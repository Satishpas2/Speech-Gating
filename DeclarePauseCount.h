//
//  DeclarePauseCount.h
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include <iostream>
#include <stdio.h>
using namespace std;
class DeclarePauseCount
{
public:
    //Power in block m
    double *powerOfEachBlock ;
    
    //Constructor
    DeclarePauseCount(int blockNumber);
    //Destructor
    ~DeclarePauseCount();
};
