//
// Created by david on 28/04/2022.
//

#include <cstring>
#include "eye.h"


void eye(int n, double m[][MAX])
{
    memset(m, 0, n*100*sizeof(double));

    for(int i = 0; i < n; i++)
        m[i][i] = 1;

}

/*
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            if( j == i)
               m[j][i] = 1;
            else
                m[j][i] = 0;
        }
    }
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            cout << m[j][i] << "\t";
        }
        cout<<endl;
    }
*/
