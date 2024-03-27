#include <bits/stdc++.h>

#include "MyMap.hpp"

using namespace std;
const int N = 210;
const int dx[] = {0, 0, -1, 1}, dy[] = {1, -1, 0, 0};

class Berth {
 private:
  int berth_id;
  int x, y;
  int transport_time;
  int loading_speed;
  stack<int> goods;
  int sum_planned, sum_gained;
  int tru_tanspost_time;
  int scheduled_end = 0;
  MyMap& mp = MyMap::GetMap();
  int direction[4] = {0, 1, 2, 3};

  struct node {
    int x, y, dis, lst_dir;
  };
  ~Berth() {
    while (!goods.empty()) goods.pop();
  }

 public:
  Berth(int id, int x, int y, int transport_time, int loading_speed) {
    this->berth_id = id;
    this->x = x, this->y = y;
    this->transport_time = transport_time;
    this->loading_speed = loading_speed;
  }
  int get_transport_time() { return transport_time; }
  float calc_berth_val(int boat_capacity) {
    return k2 * transport_time / boat_capacity +
           k3 * (goods.size() / boat_capacity) * boat_capacity / loading_speed;
  }
  void berth_bfs() {
    queue<node> q;
    bool vis[N][N];
    memset(vis, 0, sizeof(vis));
    for (int i = x; i <= x + 3; ++i)
      for (int j = y; j <= y + 3; ++j) {
        q.push({i, j, 0, -1});
        vis[i][j] = 1;
        mp.update_to_berth(berth_id, i, j, 0, -1);
        mp.update_berth(i, j, berth_id);
      }
    while (!q.empty()) {
      node now = q.front();
      q.pop();
      int pp = rand() % 4, qq = rand() % 4;
      swap(direction[pp], direction[qq]);
      for (int t = 0; t < 4; ++t) {
        int i = direction[t];
        node nxt = {now.x + dx[i], now.y + dy[i], now.dis + 1, i};
        if (mp.check_range(nxt.x, nxt.y) && mp.check_map(nxt.x, nxt.y) &&
            !vis[nxt.x][nxt.y]) {
          q.push(nxt);
          mp.update_to_berth(berth_id, nxt.x, nxt.y, nxt.dis, i ^ 1);
          vis[nxt.x][nxt.y] = 1;
        }
      }
    }
  }
};