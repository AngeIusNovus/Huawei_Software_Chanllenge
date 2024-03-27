#include <bits/stdc++.h>
using namespace std;
const int N = 210;

class MyMap {
 private:
  int n;
  char ch[N][N];
  int gds[N][N];
  int berth_map[N][N];
  int least_valued_time[N][N];
  float least_cost[N][N];
  static MyMap* mp;

  struct To_Berth {
    int dis, lst_dir;
    To_Berth() { dis = -1; }
    To_Berth(int dis, int dir) {
      this->dis = dis;
      this->lst_dir = dir;
    }
  } to_berth[10][N][N];

  void input() {
    this->n = 200;
    for (int i = 0; i < n; ++i) scanf("%s", ch[i]);
  }
  MyMap() {
    memset(gds, 0, sizeof(gds));
    memset(berth_map, -1, sizeof(berth_map));
    mp->input();
  }

 public:
  static MyMap& GetMap() {
    if (mp == NULL) mp = new MyMap();
    return *mp;
  }

  void input() {
    for (int i = 0; i < n; i++) scanf("%s", ch[i]);
  }

  int get_to_berth_dis(int id, int x, int y) { return to_berth[id][x][y].dis; }
  int get_to_berth_dir(int id, int x, int y) {
    return to_berth[id][x][y].lst_dir;
  }

  void update_berth(int x, int y, int val) { berth_map[x][y] = val; }
  void update_to_berth(int id, int x, int y, int dis, int lst_dir) {
    to_berth[id][x][y] = {dis, lst_dir};
  }
  bool check_to_berth(int id, int x, int y) {
    return to_berth[id][x][y].dis != -1;
  }
  void modify_map(int x, int y, char c) { ch[x][y] = c; }

  void modify_least_valued_time(int x, int y, int val) {
    least_valued_time[x][y] = val;
  }
  int get_least_valued_time(int x, int y) { return least_valued_time[x][y]; }
  void modify_least_cost(int x, int y, float val) { least_cost[x][y] = val; }
  float get_least_cost(int x, int y) { return least_cost[x][y]; }

  void update_goods(int x, int y, int val) { gds[x][y] = val; }

  int get_goods_id(int x, int y) { return gds[x][y]; }

  bool check_range(int x, int y) { return x >= 0 && x < n && y >= 0 && y < n; }

  bool check_map(int x, int y) { return ch[x][y] != '*' && ch[x][y] != '#'; }
};