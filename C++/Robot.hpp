#include <bits/stdc++.h>

#include "Goods_Table.hpp"
#include "MyMap.hpp"

using namespace std;

class Robot {
 private:
  int robot_id, x, y, free, sts;
  MyMap& mp = MyMap::GetMap();
  Goods_Table& GT = Goods_Table::GetGT();
  int ideal_direction;
  int real_direction;
  int has_moved;
  int dest_good;
  int self_val;
  int direction[4] = {0, 1, 2, 3};
  bool vis[N][N], p = 0;
  struct node {
    int x, y, dis, lst_dir;
  };
  struct Best_Good {
    int x, y;
    float val;
    bool operator<(const Best_Good& rhs) const& { return val > rhs.val; }
    int calc_id() { return x * 200 + y; }
  };
  priority_queue<Best_Good> best_goods;
  struct To_Robot {
    int dis, dir;
    To_Robot() { dis = 0; }
    To_Robot(int dis, int dir) {
      this->dis = dis;
      this->dir = dir;
    }
  } to_robot[N][N];

 public:
  Robot(int id) { this->robot_id = id; }
  void update(int x, int y, int free, int sts) {
    this->x = x, this->y = y;
    this->free = free;
    this->sts = sts;
  }
  void robot_bfs() {
    queue<node> q;
    p ^= 1;
    q.push({x, y, 0, -1});
    vis[x][y] = p;
    while (!q.empty()) {
      node now = q.front();
      q.pop();
      if (now.dis > 990) continue;
      if (best_goods.size() == robot_num &&
          best_goods.top().val * (mp.get_least_cost(now.x, now.y) + now.dis) >=
              200)
        continue;
      swap(direction[rand() & 3], direction[rand() & 3]);
      for (int t = 0; t < 4; ++t) {
        int i = direction[t];
        node nxt = {now.x + dx[i], now.y + dy[i], now.dis + 1, now.lst_dir};
        if (now.lst_dir == -1) nxt.lst_dir = i;
        if (mp.check_range(nxt.x, nxt.y) && mp.check_map(nxt.x, nxt.y) &&
            (vis[nxt.x][nxt.y] ^ p)) {
          q.push(nxt);
          to_robot[nxt.x][nxt.y] = {nxt.lst_dir, nxt.dis};
          vis[nxt.x][nxt.y] = p;
          if (mp.get_goods_id(nxt.x, nxt.y)) {
            Goods* good = GT.get_goods(nxt.x, nxt.y);
            if (!good->check_alive(frame_time + nxt.dis + 10)) continue;
            float ans = good->get_val(nxt.dis);
            Best_Good best_good = {nxt.x, nxt.y, ans};
            best_goods.push(best_good);
            if (best_goods.size() > robot_num) best_goods.pop();
          }
        }
      }
    }
  }
  int get_free() { return free; }
  void set_ideal_direction(int ideal_direction) {
    this->ideal_direction = ideal_direction;
  }
  void set_dest_goods(int id) { dest_good = id; }
  PFI get_best_goods() {
    while (!best_goods.empty()) {
      Best_Good goods = best_goods.top();
      best_goods.pop();
      PFI tmp = {goods.val, goods.calc_id()};
      return tmp;
    }
    return {0, -1};
  }
  int get_x() { return x; }
  int get_y() { return y; }
  int get_to_robot_dir(int x, int y) { return to_robot[x][y].dir; }
};