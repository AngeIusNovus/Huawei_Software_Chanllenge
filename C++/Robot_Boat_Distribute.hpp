#include <bits/stdc++.h>

#include "Berth.hpp"
#include "Goods_Table.hpp"
#include "MyMap.hpp"
#include "Robot.hpp"

#define PFI pair<float, int>

using namespace std;

extern int calc_id(int x, int y) { return x * 200 + y; }
const float eps = 1e-5;

class MCMF {
 private:
  struct Edge {
    int vet, nxt, f;
    double w;
  } E[25050];
  int head[40050], cur[40050], tot;
  int s = 40010;
  int t = 40011;
  double dis[40050];
  bool vis[40050];
  int tot_nodes = 0;
  int idx_x[205], idx_y[205];
  bool sim[N][N];

  bool spfa() {
    for (int i = 0; i < tot_nodes; ++i) {
      int t = calc_id(idx_x[i], idx_y[i]);
      dis[t] = 1e18;
      cur[t] = head[t];
    }
    for (int i = 40000; i <= t; ++i) {
      dis[i] = 1e18;
      cur[i] = head[i];
    }
    queue<int> q;
    double INF = 1e18;
    q.push(s), dis[s] = 0, vis[s] = 1;
    while (!q.empty()) {
      int u = q.front();
      q.pop(), vis[u] = 0;
      for (int i = head[u]; i; i = E[i].nxt) {
        int v = E[i].vet;
        if (E[i].f && dis[v] - dis[u] - E[i].w > eps) {
          dis[v] = dis[u] + E[i].w;
          if (!vis[v]) q.push(v), vis[v] = 1;
        }
      }
    }
    return dis[t] != INF;
  }

  int dfs(int u, int t, int flow) {
    if (u == t) return flow;
    vis[u] = 1;
    int res = 0;
    for (int& i = cur[u]; i && res < flow; i = E[i].nxt) {
      int v = E[i].vet;
      if (!vis[v] && E[i].f && fabs(dis[v] - dis[u] - E[i].w) < eps) {
        int x = dfs(v, t, min(E[i].f, flow - res));
        if (x) E[i].f -= x, E[i ^ 1].f += x, res += x;
      }
    }
    vis[u] = 0;
    return res;
  }

 public:
  MCMF() {
    memset(sim, 0, sizeof(sim));
    memset(head, 0, sizeof(head));
  }
  void clear() {
    tot = 1;
    for (int i = 0; i < tot_nodes; ++i) {
      sim[idx_x[i]][idx_y[i]] = 0;
      head[calc_id(idx_x[i], idx_y[i])] = 0;
    }
    for (int i = 40000; i <= t; ++i) head[i] = 0;
    tot_nodes = 0;
  }
  bool check_simed(int x, int y) { return sim[x][y]; }
  void colorize(int x, int y) {
    idx_x[tot_nodes] = x;
    idx_y[tot_nodes] = y;
    sim[x][y] = 1;
    tot_nodes++;
  }
  void add(int u, int v, int f, double val) {
    E[++tot] = {v, head[u], f, val}, head[u] = tot;
    E[++tot] = {u, head[v], 0, -val}, head[v] = tot;
  }
  void adds(int v, int f, double val) { add(s, v, f, val); }
  void addt(int u, int f, double val) { add(u, t, f, val); }

  void dinic() {
    int max_flow = 0;
    while (spfa()) {
      int x;
      while ((x = dfs(s, t, 1e9))) max_flow += x;
    }
  }
  int get_target(int u) {
    for (int i = head[u]; i; i = E[i].nxt) {
      if (!E[i].f) return E[i].vet;
    }
  }
};

class Robot_Boat_Distribute {
 private:
  int n, robot_num, berth_num, boat_num;
  int boat_capacity;
  MyMap& mp = MyMap::GetMap();
  Goods_Table& GT = Goods_Table::GetGT();
  unique_ptr<Berth> berth[10];
  unique_ptr<Robot> robot[10];
  MCMF mcmf;

  PFI get_best_berth(int x, int y) {
    float ans = 1e18;
    int dest_berth = -1;
    for (int i = 0; i < berth_num; ++i) {
      if (mp.check_to_berth(i, x, y)) continue;
      float tmp = k1 * mp.get_to_berth_dis(i, x, y) +
                  berth[i]->calc_berth_val(boat_capacity);
      if (tmp < ans) {
        ans = tmp;
        dest_berth = i;
      }
    }
    return {ans, dest_berth};
  }

 public:
  Robot_Boat_Distribute(int n, int robot_num, int berth_num, int boat_num) {
    this->n = n;
    this->robot_num = robot_num;
    this->berth_num = berth_num;
    this->boat_num = boat_num;
  }

  void Init() {
    mp = MyMap::GetMap();
    for (int i = 0; i < berth_num; i++) {
      int id, x, y, transport_time, loading_speed;
      scanf("%d", &id);
      scanf("%d%d%d%d", &x, &y, &transport_time, &loading_speed);
      berth[id] = make_unique<Berth>(
          new Berth(id, x, y, transport_time, loading_speed));
      berth[id]->berth_bfs();
    }
    for (int i = 0; i < robot_num; ++i)
      robot[i] = make_unique<Robot>(new Robot(i));
    INIT_berth();
    scanf("%d", &boat_capacity);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        if (mp.check_map(i, j)) {
          int res = 1e7;
          float cost = 1e18;
          bool flag = 1;
          for (int k = 0; k < berth_num; ++k) {
            if (mp.check_to_berth(k, i, j)) {
              flag = 0;
              cost = min(cost, (float)1.0 * berth[k]->get_transport_time() /
                                       boat_capacity +
                                   mp.get_to_berth_dis(k, i, j));
              res = min(res, berth[k]->get_transport_time() +
                                 mp.get_to_berth_dis(k, i, j));
            }
          }
          mp.modify_least_valued_time(i, j, res);
          mp.modify_least_cost(i, j, cost);
          if (flag) mp.modify_map(i, j, '#');
        }
    for (int i = 0; i < boat_num; i++) {
      boat[i].reset(i);
    }
    char okk[100];
    scanf("%s", okk);
    printf("OK\n");
    fflush(stdout);
  }

  void robot_manage() {
    for (int i = 0; i < robot_num; ++i) robot[i]->set_ideal_direction(-1);
    for (int u = 0; u < robot_num; ++u) {
      if (robot[u]->get_free()) {
        int x = robot[u]->get_x(), y = robot[u]->get_y();
        PFI tmp = get_best_berth(x, y);
        assert(tmp.second >= 0 && tmp.second < berth_num);
        robot[u]->set_ideal_direction(mp.get_to_berth_dir(tmp.second, x, y));
        continue;
      }
      int v = mcmf.get_target();
      int x = v / 200, y = v % 200;
      robot[u]->set_ideal_direction(robot[u]->get_to_robot_dir(x, y));
      robot[u]->set_dest_goods(mp.get_goods_id(x, y));
    }
  }

  int frame_input() {
    int num;
    scanf("%d", &num);
    for (int i = 1; i <= num; i++) {
      int x, y, val;
      scanf("%d%d%d", &x, &y, &val);
      if (mp.check_map(x, y)) continue;
      GT.insert(x, y, val, frame_time);
    }
    GT.Clear_Disappeared_Goods();
    mcmf.clear();
    int tot_notfree_robot = 0;
    for (int i = 0; i < robot_num; i++) {
      int x, y, free, sts;
      scanf("%d%d%d%d", &free, &x, &y, &sts);
      if (mp.get_least_valued_time(x, y) + frame_time >= dividing_frame)
        free = 1;
      robot[i]->update(x, y, free, sts);
      if (!free && sts) {
        robot[i]->robot_bfs();
        mcmf.adds(40000 + i, 1, 0);
      }
    }
    for (int u = 0; u < robot_num; ++u) {
      if (robot[u]->get_free()) continue;
      PFI good = robot[u]->get_best_goods();
      while (good.second != -1) {
        int x = good.second / 200, y = good.second % 200;
        if (!mcmf.check_simed(x, y)) {
          mcmf.colorize(x, y);
          mcmf.addt(good.second, 1, 0);
        }
        mcmf.add(40000 + u, good.second, 1, -good.first);
      }
    }
    mcmf.dinic();
    robot_manage();
    for (int i = 0; i < 5; i++) scanf("%d%d\n", &boat[i].status, &boat[i].pos);
    char okk[100];
    scanf("%s", okk);
    return frame_time;
  }
};