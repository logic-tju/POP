#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

int main(int argc, char* argv[]){
    ifstream f(argv[1]);
    map<pair<unsigned, unsigned>, unsigned> edgeMap;
    string line;
    char c;
    unsigned s, t, w;
    unsigned count = 0;
    ofstream out(argv[2]);
    while(getline(f, line)){
        istringstream is(line);
        is >> c;
        if(c != 'a') continue;
        is >> s >> t >> w;
        pair<unsigned, unsigned> temp{s, t};
        if(edgeMap.find(temp) != edgeMap.end()){
        	++count;
            //cout << "Edge: " << s << "-" << t << ":" << w << "/" << edgeMap[temp] << endl;
        }else{
            edgeMap[temp] = w;
        }
    }

	out << edgeMap.size() << endl;
	out << "Original Edges: " << edgeMap.size() + count << ". Duplicated edges: " << count << endl;
    for(auto e : edgeMap){
        out << "a " << e.first.first << ' ' << e.first.second << ' ' << e.second << endl;
    }
    

    return 0;
}
