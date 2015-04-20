#include <iostream>
#include "Option.h"
#include <queue>
using namespace std;

int main(){
    Option O(1.0,0.08,0.03,2,0.1,0.2,0.5,2.0);
    cout << O.getPrice() << endl;

    EuropeanCall E(1.0,0.08,0.03,2,0.1,0.2,0.5,2.0,1.0);
    cout << E.getPrice() << endl;
    

    return 0;
}