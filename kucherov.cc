#include <iostream>

using namespace std;


int L[200];
int U[200];
int B[200];
int parts=3;
int sigma=2;
int mm=2;



int nodesld(int l, int d, int (&nodes)[200][10])
{
    int ret=0;
 //   cout << " params " << l << " " << d << endl;
 
    if( d < 0 )
        return 0;
    
    if(l==0 && d==0){
        nodes[l][d+1] = 1;
        return 1;
    }
 
    if(d==0 && l==0){
        nodes[l][d+1] = 1;
        return 1;
    }


    // recursive case
    
    if(l >= 1 && d >= L[l] && d <= U[l]){

        int v1,v2;
        v1 = nodes[l-1][d+1];
        v2 = nodes[l-1][d];
   
        // distinguish recursive call from  lookup
        if( v1 == -1){
            nodes[l-1][d+1] = nodesld(l-1,d,nodes);
            v1 = nodes[l-1][d+1];
        }
        
        // distinguish recursive call from  lookup
        if( v2 == -1){
                nodes[l-1][d] = nodesld(l-1,d-1,nodes);
                v2 = nodes[l-1][d];
        }
        
  //      cout << l << "," << d << " addition " << v1 << " " << v2 << endl;
  
        ret = v1+(sigma-1)*v2;
        nodes[l][d+1] = ret;
        }
        else{
            nodes[l][d+1]=0;
            ret=0;
        }
    
    return ret;
}


void kucherov(int (&nl)[200], int (&nodes)[200][10]){

    
//    for (int i=1; i<=B[parts]; i++) {
//        cout << i << " : " << L[i] << " : " << U[i] << endl;
//    }
    
    
    for(int l=0; l<=B[parts]; l++){
        for(int d=0; d<=mm+1; d++)
            nodes[l][d]=-1;
    }
    
    for(int l=0; l<=B[parts]; l++){
        nl[l]=nodesld(l,L[l],nodes);
        for(int d=L[l]+1; d<=U[l]; d++)
            nl[l] += nodesld(l,d,nodes);
    }
    
}


int main(){
    
    int nl[200];
    int nodes[200][10];
    
	std::cout << "Alphabet size " << endl;
    cin >> sigma;
    
    std::cout << "Mismatches " << endl;
    cin >> mm;
    
    std::cout << "Parts " << endl;
    cin >> parts;
    
    
    
    B[0]=0;
    for (int i=1; i<=parts; i++) {
        std::cout << "Length part " << i << endl;
        cin >> B[i];
    }
    for (int i=1; i<=parts; i++) {
        B[i]=B[i]+B[i-1];
    }
    
    L[0]=U[0]=0;
    for (int i=1; i<=parts; i++) {
        std::cout << "Lower bounds part " << i << endl;
        int l;
        cin >> l;
        for(int j=B[i];j>B[i-1];j--)
            L[j]=l;
        
        std::cout << "Upper bounds part " << i << endl;
        int u;
        cin >> u;
        for(int j=B[i];j>B[i-1];j--)
            U[j]=u;
        
    }
    
    
    
    
    //
    //    for (int i=1; i<=B[parts]; i++) {
    //        cout << i << " : " << L[i] << " : " << U[i] << endl;
    //    }
    //
    // adjust bounds
    for(int i=parts; i>=1;i--){
        int j=B[i]; // upper index of part
        int l=L[j]; // the lower bound that has to be reached
        
        while( j-1 > B[i-1] && l>0 ){
            L[j-1] = l-1;
            j--;
            if( l>1 ) l--;
        }
    }
    
    kucherov(nl,nodes);


//    for(int l=0; l<=B[parts]; l++){
//        for(int d=0; d<=mm; d++)
//            cout << "Nodes[" << l << "][" << d << "]= " << nodes[l][d+1] <<  endl;
//    }
    
    int sum=0;
    for(int l=0; l<=B[parts]; l++)
        sum += nl[l];

    cout << "Nodes: " << sum << endl;

}
