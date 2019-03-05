#include <iostream>
#include <vector>
#include <string>

using namespace std;


void compute_edit_distance(vector<vector<int> >& F,string& S1, string& S2){
    for (int i=1; i <= S1.length(); i++){
        for (int j=1; j<= S2.length(); j++) {
            int d=F[i-1][j-1];
            int l=F[i-1][j];
            int u=F[i][j-1];
            
            if( S1[i-1] == S2[j-1])
                F[i][j] = d;
            else
                F[i][j] = d+1;
            if( u < F[i][j] )
                F[i][j] = u+1;
            if( l < F[i][j] )
                F[i][j] = l+1;
        }
    }
}

int main(){
    
    string S1="Hello wie geht es Dir?";
    string S2="Hallo wie gehts Dir?";
    
    vector< vector<int> > F;
    F.resize(S1.length()+1);
    for (int i=0; i <= S1.length(); i++) {
        F[i].resize(S2.length()+1);
    }
    
    for (int i=0; i<= S1.length(); i++) {
        F[i][0]=i;
    }
    for (int j=0; j<= S2.length(); j++) {
        F[0][j]=j;
    }
    
    compute_edit_distance(F,S1,S2);
    
    for (int i=0; i <= S1.length(); i++){
        cout << "i = " << i << endl;
        for (int j=0; j<= S2.length(); j++) {
            cout << F[i][j] << " ";
        }
        cout << endl;
    }
}
