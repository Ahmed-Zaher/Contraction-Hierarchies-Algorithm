typedef long long ll;
typedef pair<ll, ll> Node;
typedef vector<vector<Node>> Graph;
typedef priority_queue<Node, vector<Node>, greater<Node>> DijkstraQueue;

const ll INF = 1e15 + 10;

#define pb push_back
#define mp make_pair
#define fs first
#define sc second

class CH {
private:
	int n, m, s, t, order;
	vector<Graph> G, G_CH;
	vector<vector<ll>> dist;
	vector<DijkstraQueue> pq;
	DijkstraQueue imp_pq;
	vector<bool> contracted;
	vector<int> imp;
	vector<int> level;
	vector<int> contr_neighbours;
	void read() {
		graph >> n >> m;
		G.assign(2, Graph(n + 1, vector<Node>()));
		int u, v, w;
		for (int i = 0; i < m; ++i) {
			graph >> u >> v >> w;
			if (u == v)
				continue;
			Connect(G[0][u], v, w);
			Connect(G[1][v], u, w);
		}
	}
	void preprocess() {
		SetOrder();
		BuildG_CH();
	}
	void Connect(vector<Node>&E, int v, ll w) {
		for (auto& p : E)
			if (p.fs == v) {
				p.sc = min(p.sc, w);
				return;
			}
		E.pb(mp(v, w));
	}
	void SetOrder() {
		contracted.assign(n + 1, 0);
		imp.assign(n + 1, 0);
		level.assign(n + 1, 0);
		contr_neighbours.assign(n + 1, 0);
		for (int i = 1; i <= n; ++i)
			imp_pq.push(mp(-n, i));
		int curr_node, new_imp;
		order = 1;
		while (!imp_pq.empty()) {
			curr_node = imp_pq.top().sc;
			imp_pq.pop();
			new_imp = GetImportance(curr_node);
			if (imp_pq.empty() || new_imp - imp_pq.top().fs <= 10) {
				imp[curr_node] = order++;
				contracted[curr_node] = 1;
				ContractNode(curr_node);
			} else
				imp_pq.push(mp(new_imp, curr_node));
		}
	}
	int GetImportance(int x) {
		int u, v, shortcuts = 0, in_out = 0;
		for (int i = 0; i < 2; ++i)
			for (auto&p : G[i][x])
				if (!contracted[p.fs])
					++in_out;
		for (auto& p1 : G[1][x])
			for (auto&p2 : G[0][x]) {
				u = p1.fs;
				v = p2.fs;
				if (!contracted[u] && !contracted[v])
					++shortcuts;
			}
		int edge_diff = shortcuts - in_out;
		return edge_diff + 2 * contr_neighbours[x] + level[x];
	}
	void ContractNode(int x) {
		int u;
		ll w, mx = GetMaxEdge(x);
		set<pair<int, ll>> out_edges;
		for (auto& p : G[0][x])
			if (!contracted[p.fs])
				out_edges.insert(mp(p.fs, p.sc));
		for (auto& p : G[1][x]) {
			u = p.fs;
			if (contracted[u])
				continue;
			w = p.sc;
			Check_Witness(u, x, w, mx, out_edges, 0);
		}
		for (int i = 0; i < 2; ++i)
			for (auto& p : G[i][x]) {
				++contr_neighbours[p.fs];
				level[p.fs] = max(level[p.fs], level[x] + 1);
			}
	}
	ll GetMaxEdge(int x) {
		ll ret = 0;
		for (auto& p1 : G[1][x])
			for (auto&p2 : G[0][x])
				if (p1.fs != p2.fs && !contracted[p1.fs] && !contracted[p2.fs])
					ret = max(ret, p1.sc + p2.sc);
		return ret;
	}
	int Check_Witness(int u, int x, ll w, ll mx, set<pair<int, ll>>&out_edges,
			bool type) {
		int a, b;
		ll curr_dist, new_dist;
		DijkstraQueue D_pq;
		unordered_map<int, ll> D_dist;
		D_pq.push(mp(0, u));
		D_dist[u] = 0;
		int iter = 250 * (n - order) / n;
		while (!D_pq.empty() && iter--) {
			curr_dist = D_pq.top().fs;
			a = D_pq.top().sc;
			D_pq.pop();
			if (curr_dist > D_dist[a])
				continue;
			for (auto& p : G[0][a]) {
				new_dist = p.sc + curr_dist;
				b = p.fs;
				if (b != x && !contracted[b])
					if (D_dist.find(b) == D_dist.end() || D_dist[b] > new_dist)
						if (D_dist[b] < mx)
							D_dist[b] = new_dist, D_pq.push(mp(new_dist, b));
			}
		}
		int v, ret = 0;
		ll new_w;
		for (auto& p : out_edges) {
			v = p.fs, new_w = w + p.sc;
			if (D_dist.find(v) == D_dist.end() || D_dist.find(v)->sc > new_w) {
				++ret;
				if (!type && u != v) {
					Connect(G[0][u], v, new_w);
					Connect(G[1][v], u, new_w);
				}
			}
		}
		return ret;
	}
	void BuildG_CH() {
		G_CH.assign(2, Graph(n + 1, vector<Node>()));
		int v;
		ll w;
		for (int u = 1; u <= n; ++u)
			for (auto& p : G[0][u]) {
				v = p.fs;
				w = p.sc;
				if (imp[v] > imp[u])
					G_CH[0][u].pb(mp(v, w));
				else
					G_CH[1][v].pb(mp(u, w));
			}
	}
	ll GetDistance() {
		dist[0][s] = dist[1][t] = 0;
		ll SP = INF;
		pq[0].push(mp(0, s));
		pq[1].push(mp(0, t));
		Node front;
		int curr_node;
		ll curr_dist;
		while (!pq[0].empty() || !pq[1].empty()) {
			if (!pq[0].empty()) {
				front = pq[0].top();
				pq[0].pop();
				curr_node = front.sc;
				curr_dist = front.fs;
				if (SP >= curr_dist)
					RelaxNodeEdges(curr_node, 0);
				SP = min(SP, dist[0][curr_node] + dist[1][curr_node]);
			}
			if (!pq[1].empty()) {
				front = pq[1].top();
				pq[1].pop();
				curr_node = front.sc;
				curr_dist = front.fs;
				if (SP >= curr_dist)
					RelaxNodeEdges(curr_node, 1);
				SP = min(SP, dist[0][curr_node] + dist[1][curr_node]);
			}
		}
		if (SP == INF)
			return -1;
		return SP;
	}
	void RelaxNodeEdges(int u, int g) {
		int v;
		ll w;
		for (auto& p : G_CH[g][u]) {
			v = p.fs, w = p.sc;
			if (dist[g][v] > dist[g][u] + w) {
				dist[g][v] = dist[g][u] + w;
				pq[g].push(mp(dist[g][v], v));
			}
		}
	}
public:
	CH() {
		cout << "Preprocessing...\n";
		read();
		preprocess();
		cout << "Preprocessing done!\n";
	}
	ll Query(int _s, int _t) {
		s = _s, t = _t;
		dist.assign(2, vector<ll>(n + 1, INF));
		pq.assign(2, DijkstraQueue());
		return GetDistance();
	}
};
