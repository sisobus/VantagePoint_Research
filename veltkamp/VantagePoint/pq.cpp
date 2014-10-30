#include <cstdio>
#include <queue>
using namespace std;

int main() {
    priority_queue<int,vector<int>,greater<int> > pq;
    for ( int i = 0 ; i < 10 ; i++ ) 
        pq.push(i);
    while ( !pq.empty() ) {
        int now = pq.top();pq.pop();
        printf("%d ",now);
    }
    return 0;
}
