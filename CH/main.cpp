#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <fstream>
#include <ctime>
#include <sstream>

using namespace std;



ifstream graph("NY_graph.txt");
ifstream queries("Queries.txt");

#include "ch.h"

int main() {
	ios::sync_with_stdio(0);
	int T = clock();
	CH G;
	cout << "Preprocessing time: " << 1.0 * (clock() - T) / CLOCKS_PER_SEC << "s\n";
	int q, s, t;
	double tot_time = 0, query_time;
	queries >> q;
	for (int i = 0; i < q; ++i) {
		queries >> s >> t;
		T = clock();
		cout << "Shortest path distance: " << G.Query(s, t) << '\n';
		query_time = 1.0 * (clock() - T) / CLOCKS_PER_SEC;
		cout << "Time: " << query_time << "s\n";
		tot_time += query_time;
	}
	cout << "Average query time: " << tot_time / q << "s\n";
	cout << "Press Enter to terminate...";
	cin.ignore();
	string term;
	getline(cin,term);
	return 0;
}
