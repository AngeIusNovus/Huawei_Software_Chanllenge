#include <bits/stdc++.h>
#pragma GCC optimize(3, "Ofast", "inline")
#define PDI pair<double, int>
#define double float
using namespace std;

const double k1 = 1.0, k2 = 1.50, k3 = 1.0;
double k4 = 2.0;
double stp = 2.0;
double MB = 4.0;
double cA, cB;

double getk4(int x) { return cB / (x + cA); }

const int n = 200;
const int robot_num = 10;
const int berth_num = 10;
const int boat_num = 5;
const int N = 210;
const int dividing_frame = 14996;
const double eps = 1e-5;

const int dx[] = {0, 0, -1, 1}, dy[] = {1, -1, 0, 0};
int direction[] = {0, 1, 2, 3};

FILE *filepointer;
int money, boat_capacity, frame_time, ideal_money, planned_money, all_money;
char ch[N][N];
int gds[N][N], berth_map[N][N];
double least_cost[N][N];
int least_valued_time[N][N];
int tot_goods = 0;
void print_berth_situation();

struct node {
  int x, y, dis, lst_dir;
};

struct To_Robot {
  int nxt_dir, dis;
  To_Robot() { nxt_dir = -1; }
  To_Robot(int dir, int ds) { nxt_dir = dir, dis = ds; }
} to_robot[robot_num][N][N];

struct To_Berth {
  int nxt_dir, dis;
  To_Berth() { dis = -1; }
  To_Berth(int b, int c) {
    nxt_dir = b;
    dis = c;
  }
  bool operator<(const To_Berth &rhs) const & { return dis < rhs.dis; }
} to_berth[10][N][N];

int get_berth(int x, int y) { return berth_map[x][y]; }

bool check_range(int x, int y) { return x >= 0 && x < n && y >= 0 && y < n; }

bool check_map(int x, int y) { return ch[x][y] != '*' && ch[x][y] != '#'; }

int calc_id(int x, int y) { return x * n + y; }

struct Goods {
  int goods_id, x, y, val, gen_time, tag, dest_berth;
  double cost;
  Goods(int id, int xx, int yy, int vval, int ttag, int dest, double c) {
    goods_id = id;
    x = xx, y = yy, val = vval;
    gen_time = frame_time;
    gds[x][y] = goods_id;
    tag = ttag;
    dest_berth = dest;
    cost = c;
  }
  void clear() {
    gds[x][y] = 0;
    tag = 1;
  }
  bool check_alive(int x) { return x - gen_time <= 1000; }
};

deque<Goods> goods;

struct Best_Good {
  int x, y;
  double val;
  bool operator<(const Best_Good &rhs) const & { return val > rhs.val; }
};

priority_queue<Best_Good> best_goods[robot_num];

struct Berth {
  int berth_id;
  int x, y;
  int transport_time;
  int loading_speed;
  stack<int> goods;
  int sum_planned;
  int sum_gained;
  int tru_tanspost_time;
  int scheduled_end = 0;
  void berth_bfs() {
    queue<node> q;
    bool vis[N][N];
    memset(vis, 0, sizeof(vis));
    for (int i = x; i <= x + 3; ++i)
      for (int j = y; j <= y + 3; ++j) {
        q.push({i, j, 0, -1});
        vis[i][j] = 1;
        to_berth[berth_id][i][j] = {-1, 0};
        berth_map[i][j] = berth_id;
      }
    while (!q.empty()) {
      node now = q.front();
      q.pop();
      int pp = rand() % 4, qq = rand() % 4;
      swap(direction[pp], direction[qq]);
      for (int t = 0; t < 4; ++t) {
        int i = direction[t];
        node nxt = {now.x + dx[i], now.y + dy[i], now.dis + 1, i};
        if (check_range(nxt.x, nxt.y) && check_map(nxt.x, nxt.y) &&
            !vis[nxt.x][nxt.y]) {
          q.push(nxt);
          to_berth[berth_id][nxt.x][nxt.y] = {i ^ 1, nxt.dis};
          vis[nxt.x][nxt.y] = 1;
        }
      }
    }
  }
} berth[berth_num + 10];

struct Robot {
  int robot_id;
  int x, y, free;
  int status;
  int ideal_direction;
  int real_direction;
  int has_moved;
  int dest_good;
  int self_val;
  bool vis[N][N], p = 0;
  Robot() { memset(vis, 0, sizeof(vis)); }
  void robot_bfs() {
    queue<node> q;
    p ^= 1;
    q.push({x, y, 0, -1});
    vis[x][y] = p;
    while (!q.empty()) {
      node now = q.front();
      q.pop();
      if (now.dis > 990) continue;
      if (best_goods[robot_id].size() == robot_num &&
          best_goods[robot_id].top().val *
                  (least_cost[now.x][now.y] + now.dis) >=
              200)
        continue;
      swap(direction[rand() & 3], direction[rand() & 3]);
      for (int t = 0; t < 4; ++t) {
        int i = direction[t];
        node nxt = {now.x + dx[i], now.y + dy[i], now.dis + 1, now.lst_dir};
        if (now.lst_dir == -1) nxt.lst_dir = i;
        if (check_range(nxt.x, nxt.y) && check_map(nxt.x, nxt.y) &&
            (vis[nxt.x][nxt.y] ^ p)) {
          q.push(nxt);
          to_robot[robot_id][nxt.x][nxt.y] = {nxt.lst_dir, nxt.dis};
          vis[nxt.x][nxt.y] = p;
          if (gds[nxt.x][nxt.y]) {
            Goods good = goods[gds[nxt.x][nxt.y] - goods.front().goods_id];
            if (!good.check_alive(frame_time + nxt.dis + 10)) continue;
            double ans = 1.0 * good.val / (good.cost + k4 * nxt.dis);
            Best_Good best_good = {nxt.x, nxt.y, ans};
            best_goods[robot_id].push(best_good);
            if (best_goods[robot_id].size() > robot_num)
              best_goods[robot_id].pop();
          }
        }
      }
    }
  }
} robot[robot_num];

int get_robot_goodsval(int id) { return robot[id].self_val; }

double calc_bb_val(int dis, int transport_time, int Size, int loading_speed) {
  return k1 * dis + k2 * transport_time / boat_capacity +
         k3 * ((Size / boat_capacity) * boat_capacity + 1) / loading_speed;
}

PDI get_best_berth(int x, int y) {
  double ans = 1e18;
  int dest_berth = -1;
  for (int i = 0; i < berth_num; ++i) {
    To_Berth now = to_berth[i][x][y];
    if (now.dis == -1) continue;
    double tmp = calc_bb_val(now.dis, berth[i].transport_time,
                             berth[i].goods.size(), berth[i].loading_speed);
    if (tmp < ans) {
      ans = tmp;
      dest_berth = i;
    }
  }
  return {ans, dest_berth};
}

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
    for (int &i = cur[u]; i && res < flow; i = E[i].nxt) {
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

  void manage() {
    for (int i = 0; i < robot_num; ++i) robot[i].ideal_direction = -1;
    for (int u = 0; u < robot_num; ++u) {
      if (robot[u].free) {
        PDI tmp = get_best_berth(robot[u].x, robot[u].y);
        // cerr << "dest berth : " << tmp.second << endl;
        assert(tmp.second >= 0 && tmp.second < berth_num);
        robot[u].ideal_direction =
            to_berth[tmp.second][robot[u].x][robot[u].y].nxt_dir;
        /*filepointer = fopen("warnings.out", "a");
        fprintf(filepointer,
                "frame_time:%d robot(%d, %d) dest_berth : %d-th(%d, %d), "
                "direct : %d\n",
                frame_time, robot[u].x, robot[u].y, tmp.second,
                berth[tmp.second].x, berth[tmp.second].y,
                robot[u].ideal_direction);
        fclose(filepointer);*/
        continue;
      }
      for (int i = head[u + 40000]; i; i = E[i].nxt) {
        if (E[i].f) continue;
        int x = E[i].vet / 200, y = E[i].vet % 200;
        robot[u].ideal_direction = to_robot[u][x][y].nxt_dir;
        robot[u].dest_good = gds[x][y];
        break;
      }
    }
  }

} mcmf;

struct Boat {
  int volume, pos, status, target, id, goods_sum;
  void reset(int id) {
    this->volume = boat_capacity;
    this->pos = -1;
    this->status = 1;
    this->target = -1;
    this->id = id;
  }
  void act() {
    if (pos == -1 && status != 0) {
      volume = boat_capacity;
      planned_money += goods_sum;
      goods_sum = 0;
    }
    if (status != 0 &&
        ((pos == -1) || (berth[pos].goods.empty() || volume == 0))) {
      if (target != -1)
        printf("ship %d %d\n", id, target);
      else if (pos != -1)
        printf("go %d\n", id);
      // print_berth_situation();
      /*filepointer = fopen("ship_command.out", "a");
      if (target != -1)
        fprintf(filepointer,
                "Case3: frametime:%d ship_status:%d ship_pos:%d volume:%d "
                "target:%d COMMAND:ship %d %d\n",
                frame_time, status, pos, volume, target, id, target);
      else if (pos != -1)
        fprintf(filepointer,
                "Case4: frametime:%d ship_status:%d ship_pos:%d volume:%d "
                "target:%d COMMAND:go %d\n",
                frame_time, status, pos, volume, target, id);
      fclose(filepointer);*/
    } else if (status == 1 && pos != -1 && !berth[pos].goods.empty() &&
               volume != 0) {
      int tmp = berth[pos].loading_speed;
      while (tmp != 0 && volume != 0 && !berth[pos].goods.empty()) {
        ideal_money += berth[pos].goods.top();
        goods_sum += berth[pos].goods.top();
        tmp--;
        volume--;
        if (berth[pos].sum_planned > 0) berth[pos].sum_planned--;
        berth[pos].goods.pop();
      }
    }
  }
} boat[10];

void print_berth_situation() {
  for (int i = 0; i < berth_num; i++) {
    filepointer = fopen("berth_situation.out", "a");
    fprintf(filepointer,
            "Frame:%d Berth:%d Sum:%d Gained:%d Planned:%d scheduled_end:%d\n",
            frame_time, i, berth[i].goods.size(), berth[i].sum_gained,
            berth[i].sum_planned, berth[i].scheduled_end);
    fclose(filepointer);
  }
  filepointer = fopen("berth_situation.out", "a");
  fprintf(filepointer,
          "Frame:%d money:%d ideal_money:%d planned_money:%d all_money:%d\n",
          frame_time, money, ideal_money, planned_money, all_money);
  fclose(filepointer);
  return;
}

void INIT_berth() {
  int mn = 2e9, relay_id = -1;
  for (int i = 0; i < berth_num; i++) {
    if (berth[i].transport_time < mn) {
      mn = berth[i].transport_time;
      relay_id = i;
    }
  }
  for (int i = 0; i < berth_num; i++) {
    berth[i].tru_tanspost_time = berth[i].transport_time;
  }
  return;
}

void Init() {
  for (int i = 0; i < n; i++) scanf("%s", ch[i]);
  for (int i = 0; i < berth_num; i++) {
    int id;
    scanf("%d", &id);
    berth[id].berth_id = id;
    scanf("%d%d%d%d", &berth[id].x, &berth[id].y, &berth[id].transport_time,
          &berth[id].loading_speed);
    berth[id].berth_bfs();
  }
  INIT_berth();
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (check_map(i, j)) {
        least_cost[i][j] = 1e18;
        least_valued_time[i][j] = 1e7;
        bool flag = 1;
        for (int k = 0; k < berth_num; ++k) {
          if (to_berth[k][i][j].dis != -1) {
            flag = 0;
            least_cost[i][j] =
                min(least_cost[i][j],
                    (float)1.0 * berth[k].transport_time / boat_capacity +
                        to_berth[k][i][j].dis);
            least_valued_time[i][j] =
                min(least_valued_time[i][j],
                    berth[k].transport_time + to_berth[k][i][j].dis);
          }
        }
        if (flag) ch[i][j] = '#';
      }
  scanf("%d", &boat_capacity);
  cA = MB * boat_capacity / (stp - 1);
  cB = stp * cA;
  for (int i = 0; i < boat_num; i++) {
    boat[i].reset(i);
  }
  char okk[100];
  scanf("%s", okk);
  printf("OK\n");
  fflush(stdout);
}

int Input() {
  scanf("%d%d", &frame_time, &money);
  int num;
  scanf("%d", &num);
  for (int i = 1; i <= num; i++) {
    int x, y, val;
    scanf("%d%d%d", &x, &y, &val);
    if (ch[x][y] == '#') continue;
    PDI tmp = get_best_berth(x, y);
    // PDI tmp = {1, 0};
    goods.push_back({++tot_goods, x, y, val, 0, tmp.second, tmp.first});
  }
  while (!goods.empty() &&
         (!goods.front().check_alive(frame_time) || goods.front().tag)) {
    goods.front().clear();
    goods.pop_front();
  }
  double tot_port_goods = 0;
  for (int i = 0; i < berth_num; ++i) {
    tot_port_goods += berth[i].goods.size();
  }
  tot_port_goods = min(MB * boat_capacity, tot_port_goods);
  k4 = getk4(tot_port_goods);
  // cerr << "before good" << endl;
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d before good\n", frame_time);
  fclose(filepointer);*/
  mcmf.clear();
  //  cerr << "after good" << endl;
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d after pull\n", frame_time);
  fclose(filepointer);*/
  int tot_notfree_robot = 0;
  for (int i = 0; i < robot_num; i++) {
    int sts;
    robot[i].robot_id = i;
    scanf("%d%d%d%d", &robot[i].free, &robot[i].x, &robot[i].y, &sts);
    tot_notfree_robot += robot[i].free ^ 1;
    if (least_valued_time[robot[i].x][robot[i].y] + frame_time >=
        dividing_frame)
      robot[i].free = 1;
    if (!robot[i].free && sts) {
      robot[i].robot_bfs();
      mcmf.adds(40000 + i, 1, 0);
    }
  }
  /*for (int i = 0; i < robot_num)
     cerr << "after bfs" << endl;*/
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d after bfs\n", frame_time);
  fclose(filepointer);*/
  for (int u = 0; u < robot_num; ++u) {
    if (robot[u].free) continue;
    while (!best_goods[u].empty()) {
      auto good = best_goods[u].top();
      best_goods[u].pop();
      int x = good.x, y = good.y;
      if (!mcmf.check_simed(x, y)) {
        mcmf.colorize(x, y);
        mcmf.addt(calc_id(x, y), 1, 0);
      }
      mcmf.add(40000 + u, calc_id(x, y), 1, -good.val);
    }
  }
  // cerr << "before manage" << endl;
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d before mcmf\n", frame_time);
  fclose(filepointer);*/
  mcmf.dinic();
  mcmf.manage();
  //  cerr << "after manage" << endl;
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d after mcmf\n", frame_time);
  fclose(filepointer);*/
  for (int i = 0; i < 5; i++) scanf("%d%d\n", &boat[i].status, &boat[i].pos);
  char okk[100];
  scanf("%s", okk);
  return frame_time;
}

struct cod {  // 坐标类型
  int x, y, type, direct;
  pair<int, int> get_status(int x, int y) {
    if (ch[x][y] == '#' || ch[x][y] == '*') return {2, -1};
    for (int i = 0; i < robot_num; i++) {
      if (robot[i].x == x && robot[i].y == y) {
        return {robot[i].free, robot[i].ideal_direction};
      }
    }
    return {-1, -1};
  }
  bool can_move() {
    if (!check_range(x, y)) return false;
    if (type == 2) return false;
    for (int i = 0; i < robot_num; i++) {
      if (robot[i].x == x && robot[i].y == y) {
        if (robot[i].ideal_direction == -1) return false;
        return !robot[i].has_moved;
      }
    }
    return true;
  }
  cod(int x, int y) {
    this->x = x;
    this->y = y;
    pair<int, int> tmp = get_status(x, y);
    this->type = tmp.first;
    this->direct = tmp.second;
  }
};

void truly_move(cod A, int D) {
  for (int i = 0; i < robot_num; i++) {
    if (robot[i].x == A.x && robot[i].y == A.y) {
      robot[i].has_moved = true;
      robot[i].x += dx[D];
      robot[i].y += dy[D];
      /*filepointer = fopen("robot_command.out", "w");
      fprintf(filepointer, "frame_time:%d move %d %d\n", frame_time, i, D);
      fclose(filepointer);*/
      printf("move %d %d\n", i, D);
      robot[i].real_direction = D;
      return;
    }
  }
  return;
}

bool avoid(cod A, int D);

bool move(cod A, int D) {
  /*filepointer = fopen("move_avoid_information.out", "a");
  fprintf(filepointer,
          "frame_time:%d Case:move x:%d y:%d type:%d direct:%d D:%d\n",
          frame_time, A.x, A.y, A.type, A.direct, D);
  fclose(filepointer);*/
  if (!A.can_move()) return false;
  if (A.type == -1) return true;
  cod B(A.x + dx[D], A.y + dy[D]);
  if (!B.can_move()) return false;
  if (B.type == -1) {
    truly_move(A, D);
    return true;
  }
  if (A.type > B.type) {
    if (dx[D] + dx[B.direct] == 0 && dy[D] + dy[B.direct] == 0) {
      if (avoid(B, D)) {
        truly_move(A, D);
        return true;
      }
    } else {
      if (move(B, B.direct)) {
        truly_move(A, D);
        return true;
      } else if (avoid(B, D)) {
        truly_move(A, D);
        return true;
      }
    }
  }
  return false;
}

bool avoid_helper1(cod A, int D) {
  /*filepointer = fopen("move_avoid_information.out", "a");
  fprintf(filepointer,
          "frame_time:%d Case:avoid_helper1 x:%d y:%d type:%d direct:%d D:%d\n",
          frame_time, A.x, A.y, A.type, A.direct, D);
  fclose(filepointer);*/
  cod B(A.x + dx[D], A.y + dy[D]);
  if (B.type <= A.type && avoid(B, D)) {
    truly_move(A, D);
    return true;
  }
  return false;
}

bool avoid_helper2(cod A, int D) {
  /*filepointer = fopen("move_avoid_information.out", "a");
  fprintf(filepointer,
          "frame_time:%d Case:avoid_helper2 x:%d y:%d type:%d direct:%d D:%d\n",
          frame_time, A.x, A.y, A.type, A.direct, D);
  fclose(filepointer);*/
  cod B(A.x + dx[D], A.y + dy[D]);
  if (avoid(B, D)) {
    truly_move(A, D);
    return true;
  }
  return false;
}

bool avoid(cod A, int D) {
  /*filepointer = fopen("move_avoid_information.out", "a");
  fprintf(filepointer,
          "frame_time:%d Case:avoid x:%d y:%d type:%d direct:%d D:%d\n",
          frame_time, A.x, A.y, A.type, A.direct, D);
  fclose(filepointer);*/
  if (!A.can_move()) return false;
  if (A.type == -1) return true;
  int D1 = D ^ 2, D2 = D ^ 3;
  if (avoid_helper1(A, D1)) return true;
  if (avoid_helper1(A, D2)) return true;
  if (avoid_helper1(A, D)) return true;
  if (avoid_helper2(A, D1)) return true;
  if (avoid_helper2(A, D2)) return true;
  if (avoid_helper2(A, D)) return true;
  return false;
}

void Reset_status() {
  for (int i = 0; i < robot_num; i++) {
    robot[i].has_moved = false;
  }
  return;
}

void update_berth_situation(int berth_id, int goodsval) {
  berth[berth_id].goods.push(goodsval);
  berth[berth_id].sum_gained++;
  return;
}

void make_robot_command() {
  for (int i = 0; i < robot_num; i++) {
    if (robot[i].free) {
      move(cod(robot[i].x, robot[i].y), robot[i].ideal_direction);
      /*filepointer = fopen("warnings.out", "a");
      fprintf(filepointer,
              "frame_time:%d robot(%d, %d) "
              "direct : %d already moved\n",
              frame_time, robot[i].x, robot[i].y, robot[i].ideal_direction);
      fclose(filepointer);*/
    }
  }
  for (int i = 0; i < robot_num; i++) {
    if (!robot[i].free) {
      move(cod(robot[i].x, robot[i].y), robot[i].ideal_direction);
      /*filepointer = fopen("warnings.out", "a");
      fprintf(filepointer,
              "frame_time:%d robot(%d, %d) "
              "direct : %d already moved\n",
              frame_time, robot[i].x, robot[i].y, robot[i].ideal_direction);
      fclose(filepointer);*/
    }
  }
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d all_robots already moved\n", frame_time);
  fclose(filepointer);*/
  // cerr << "already moved" << endl;
  for (int i = 0; i < robot_num; ++i)
    if ((!robot[i].free) &&
        (gds[robot[i].x][robot[i].y] == robot[i].dest_good) &&
        (robot[i].dest_good != 0)) {
      /*filepointer = fopen("robot_command.out", "w");
      fprintf(filepointer, "frame_time:%d get %d\n", frame_time, i);
      fclose(filepointer);*/
      printf("get %d\n", i);
      auto good = goods[gds[robot[i].x][robot[i].y] - goods.front().goods_id];
      good.clear();
      robot[i].self_val = good.val;
      all_money += good.val;
    }
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d already get\n", frame_time);
  fclose(filepointer);*/
  // cerr << "already get" << endl;
  bool flag = 1;
  for (int i = 0; i < robot_num; ++i)
    if (robot[i].free && ch[robot[i].x][robot[i].y] == 'B') {
      /*filepointer = fopen("robot_command.out", "w");
      fprintf(filepointer, "frame_time:%d pull %d\n", frame_time, i);
      fclose(filepointer);*/
      printf("pull %d\n", i);
      flag = 0;
      update_berth_situation(get_berth(robot[i].x, robot[i].y),
                             get_robot_goodsval(i));
    }
  /*filepointer = fopen("warnings.out", "a");
  fprintf(filepointer, "frame_time:%d already pull\n", frame_time);
  fclose(filepointer);*/
  // cerr << "already pull" << endl;
  if (flag) return;
  for (auto good : goods) {
    if (good.tag) continue;
    PDI tmp = get_best_berth(good.x, good.y);
    good.cost = tmp.first;
    good.dest_berth = tmp.second;
  }
  return;
}

struct Boat_command {
  int boat_id;
  int target;
  double val;
  int reach_time = 1;
  int pick_time;
  int ideal_got_sum;
  int waiting_time;
  int start_berth_id = -1;

  Boat_command(int boat_id, int target) {
    this->boat_id = boat_id;
    this->target = target;
  }
  bool operator<(const Boat_command &rhs) const & { return val > rhs.val; }
  void make_true() {
    if (start_berth_id != boat[boat_id].pos)
      /*berth[target].sum_planned += ideal_got_sum,
          berth[target].scheduled_end =
              reach_time + waiting_time + pick_time - 1,*/
      target = -1;
    if (target == -1) {
      boat[boat_id].target = target;
    } else {
      boat[boat_id].target = target;
      berth[target].sum_planned += ideal_got_sum;
      int end = reach_time + waiting_time + pick_time - 1;
      berth[target].scheduled_end = end;
    }
    return;
  }
};

int dis(int start_id, int end_id) {
  if (start_id == -1) return berth[end_id].tru_tanspost_time;
  if (end_id == -1) return berth[start_id].tru_tanspost_time;
  return 500;
}

int get_waiting_time(int reach_time, int berth_id) {
  return max(0, berth[berth_id].scheduled_end - reach_time + 1);
}

const double ideal_K = 0.4;  // 船只对未来预期过好，需要给予抑制

Boat_command evaluate_boat_command(int boat_id, int start_berth_id,
                                   int berth_id, int corrected_val = 0) {
  Boat_command boat_command(boat_id, berth_id);
  double produce_rate = 1.0 * berth[berth_id].sum_gained / frame_time;
  int backing_time = berth[berth_id].tru_tanspost_time,
      sailing_time = dis(start_berth_id, berth_id);
  double all_produce_rate = 0;
  for (int i = 0; i < berth_num; i++) {
    all_produce_rate += 1.0 * berth[i].sum_gained / frame_time;
  }
  int ideal_got_sum =
      min(boat[boat_id].volume,
          (int)((sailing_time + corrected_val) * produce_rate * ideal_K -
                berth[berth_id].sum_planned + berth[berth_id].goods.size() -
                all_produce_rate * corrected_val / 5));
  int waiting_time =
      get_waiting_time(frame_time + sailing_time + corrected_val, berth_id);
  double pick_time =
      1.0 * ideal_got_sum / berth[berth_id].loading_speed + waiting_time;
  boat_command.val = 1.0 * ideal_got_sum /
                     (pick_time + sailing_time + backing_time + corrected_val);
  int really_cost_time =
      frame_time + pick_time + sailing_time + backing_time + corrected_val;
  /*if (corrected_val != 0)
    really_cost_time = frame_time + berth[boat_id].tru_tanspost_time;*/
  if (really_cost_time > dividing_frame) boat_command.val = 0;
  boat_command.pick_time = pick_time;
  boat_command.reach_time = frame_time + sailing_time;
  boat_command.ideal_got_sum = ideal_got_sum;
  boat_command.waiting_time = waiting_time;
  boat_command.start_berth_id = start_berth_id;
  /*filepointer = fopen("boat_evaluate.out", "a");
  fprintf(filepointer,
          "frame_time:%d boat_id:%d val:%.7lf fenzi:%d fenmu:%.7lf "
          "start_berth:%d end_berth:%d corrected_val:%d\n",
          frame_time, boat_id, boat_command.val, ideal_got_sum,
          (pick_time + sailing_time + backing_time + corrected_val),
          start_berth_id, berth_id, corrected_val);
  fclose(filepointer);*/
  return boat_command;
}

void make_boat_command() {
  if (frame_time == 1) {
    vector<PDI> vec;
    for (int i = 0; i < berth_num; i++) {
      vec.push_back({1.0 * berth[i].tru_tanspost_time / boat_capacity +
                         1.0 / berth[i].loading_speed,
                     i});
    }
    sort(vec.begin(), vec.end());
    for (int i = 0; i < boat_num; i++) {
      Boat_command tmp(i, vec[i].second);
      tmp.reach_time = berth[vec[i].second].tru_tanspost_time;
      tmp.make_true();
      boat[i].act();
    }
    return;
  }
  for (int i = 0; i < boat_num; i++) {
    if (boat[i].status != 0 && boat[i].pos == -1) {
      vector<Boat_command> vec;
      for (int j = 0; j < berth_num; j++) {
        vec.push_back(evaluate_boat_command(i, boat[i].pos, j, 0));
      }
      sort(vec.begin(), vec.end());
      if (vec[0].val < eps) vec[0].target = -1;
      vec[0].make_true();
    } else if (boat[i].status != 0 && boat[i].pos != -1 &&
               (berth[boat[i].pos].goods.empty() || boat[i].volume == 0)) {
      vector<Boat_command> vec;
      for (int j = 0; j < berth_num; j++) {
        vec.push_back(evaluate_boat_command(
            i, -1, j, berth[boat[i].pos].tru_tanspost_time));
        if (boat[i].pos != j)
          vec.push_back(evaluate_boat_command(i, boat[i].pos, j, 0));
      }
      sort(vec.begin(), vec.end());
      if (vec[0].val < eps) vec[0].target = -1;
      vec[0].make_true();
    }
    boat[i].act();
  }
  return;
}
int main() {
  // freopen("test.in", "r", stdin);
  // freopen("test.out", "w", stdout);
  srand(0);
  /*filepointer = fopen("ship_command.out", "w");
  fprintf(filepointer, "START!!!\n");
  fclose(filepointer);
  filepointer = fopen("berth_situation.out", "w");
  fprintf(filepointer, "START!!!\n");
  fclose(filepointer);
  filepointer = fopen("warnings.out", "w");
  fprintf(filepointer, "START!!!\n");
  fclose(filepointer);
  filepointer = fopen("robot_command.out", "w");
  fprintf(filepointer, "START!!!\n");
  fclose(filepointer);
  filepointer = fopen("bb_information.out", "w");
  fprintf(filepointer, "START!!!\n");
  fclose(filepointer);
  filepointer = fopen("boat_evaluate.out", "w");
  fprintf(filepointer, "START!!!\n");
  fclose(filepointer);
  filepointer = fopen("move_avoid_information.out", "w");
  fprintf(filepointer, "START!!!\n");
  fclose(filepointer);
  filepointer = fopen("bb_information.out", "w");
  for (int i = 0; i < berth_num; i++) {
    fprintf(filepointer,
            "berth_id:%d loading_speed:%d transport_time:%d "
            "tru_transport_time:%d \n",
            berth[i].berth_id, berth[i].loading_speed, berth[i].transport_time,
            berth[i].tru_tanspost_time);
  }
  fprintf(filepointer, "boat_capacity: %d\n", boat_capacity);
  fclose(filepointer);*/
  Init();
  for (int zhen = 1; zhen <= 15000; zhen++) {
    int id = Input();
    // puts("OK");
    // for (int i = 0; i < robot_num; ++i) robot[i].ideal_direction = rand() &
    // 3;
    /*filepointer = fopen("warnings.out", "a");
    fprintf(filepointer, "frame_time:%d already managed\n", frame_time);
    fclose(filepointer);*/
    // cerr << "already managed" << endl;
    make_robot_command();
    /*filepointer = fopen("warnings.out", "a");
    fprintf(filepointer, "frame_time:%d next step\n", frame_time);
    fclose(filepointer);*/
    // cerr << "next step" << endl;
    make_boat_command();
    /*filepointer = fopen("warnings.out", "a");
    fprintf(filepointer, "frame_time:%d boat already moved\n", frame_time);
    fclose(filepointer);*/
    // cerr << "boat already moved" << endl;
    puts("OK");
    Reset_status();
    fflush(stdout);

    if (frame_time >= 14999) {
      filepointer = fopen("warnings.out", "a");
      fprintf(filepointer, "stp:%.2f frame_time:%d zhen:%d money:%d \n", stp,
              frame_time, zhen, money);
      fclose(filepointer);
      print_berth_situation();
    }
    if (frame_time >= 14999) break;
  }

  return 0;
}
