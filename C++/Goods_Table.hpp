#include <bits/stdc++.h>

#include "MyMap.hpp"

using namespace std;

class Goods {
 private:
  int goods_id, x, y, gen_time, val, tag, dest_berth;
  float cost;

 public:
  Goods(int id, int x, int y, int val, int gen_time) {
    this->goods_id = id;
    this->x = x, this->y = y;
    this->gen_time = gen_time;
    this->val = val;
    tag = 0, dest_berth = -1;
  }
  ~Goods() {}
  void clear() {
    MyMap &mp = MyMap::GetMap();
    mp.update_goods(x, y, 0);
    tag = 1;
  }
  bool check_alive(int x) { return (x - gen_time <= 1000) && (!tag); }
  int get_goods_id() { return goods_id; }
  void update(int dest_berth, float cost) {
    this->dest_berth = dest_berth;
    this->cost = cost;
  }
  int get_val(int dis) { return 1.0 * val / (k4 * dis + cost); }
};

class Goods_Table {
 private:
  int tot_goods;
  deque<Goods *> GT;
  static Goods_Table *goods_table;
  MyMap &mp = MyMap::GetMap();
  Goods_Table() {
    tot_goods = 0;
    GT.clear();
  }

 public:
  ~Goods_Table() {
    for (auto good : GT) {
      delete (good);
    }
  }
  static Goods_Table &GetGT() {
    if (goods_table == NULL) goods_table = new Goods_Table();
    return *goods_table;
  }
  void insert(int x, int y, int val, int gen_time) {
    GT.push_back(new Goods(++tot_goods, x, y, val, gen_time));
    MyMap &mp = MyMap::GetMap();
    mp.update_goods(x, y, tot_goods);
  }
  void Clear_Disappeared_Goods() {
    while (!GT.empty() && (!GT.front()->check_alive(frame_time))) {
      GT.front()->clear();
      GT.pop_front();
    }
  }
  Goods *get_goods(int x, int y) {
    if (!mp.get_goods_id(x, y)) return NULL;
    return GT[mp.get_goods_id(x, y) - GT.front()->get_goods_id()];
  }
};