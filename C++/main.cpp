#include <bits/stdc++.h>
using namespace std;

const int n = 200;
const int robot_num = 10;
const int berth_num = 4;
const int N = 210;
const int dx[] = {0, 0, -1, 1}, dy[] = {1, -1, 0, 0};

int money, boat_capacity, frame_id;
char ch[N][N];
int gds[N][N];

struct node {
  int x, y, dis, lst_dir;
};

struct To_Robot {
  int dest_bot, nxt_dir, dis;
  bool operator<(const To_Robot &rhs) const & { return dis < rhs.dis; }
};

vector<To_Robot> to_robot[N][N];

struct To_Berth {
  int dest_berth, nxt_dir, dis;
  bool operator<(const To_Berth &rhs) const & { return dis < rhs.dis; }
};

vector<To_Berth> to_berth[N][N];

bool check_range(int x, int y) { return x >= 0 && x < n && y >= 0 && y < n; }

bool check_map(int x, int y) { return ch[x][y] != '*' && ch[x][y] != '#'; }

struct Robot {
  int robot_id;
  int x, y, goods;
  int status;
  int mbx, mby;
  Robot() {}
  Robot(int startX, int startY) {
    x = startX;
    y = startY;
  }
  void robot_bfs() {
    queue<node> q;
    bool vis[N][N];
    memset(vis, 0, sizeof(vis));
    q.push({x, y, 0, -1});
    vis[x][y] = 1;
    while (!q.empty()) {
      node now = q.front();
      q.pop();
      for (int i = 0; i < 4; ++i) {
        node nxt = {now.x + dx[i], now.y + dy[i], now.dis + 1, i};
        if (check_range(nxt.x, nxt.y) && check_map(nxt.x, nxt.y) &&
            !vis[nxt.x][nxt.y]) {
          q.push(nxt);
          to_robot[nxt.x][nxt.y].push_back({robot_id, i, nxt.dis});
          vis[nxt.x][nxt.y] = 1;
        }
      }
    }
  }
} robot[robot_num + 10];

struct Berth {
  int berth_id;
  int x, y;
  int transport_time;
  int loading_speed;
  void berth_bfs() {
    queue<node> q;
    bool vis[N][N];
    memset(vis, 0, sizeof(vis));
    for (int i = x; i <= x + 3; ++i)
      for (int j = y; j <= y + 3; ++j) q.push({i, j, 0, -1}), vis[i][j] = 1;
    while (!q.empty()) {
      node now = q.front();
      q.pop();
      for (int i = 0; i < 4; ++i) {
        node nxt = {now.x + dx[i], now.y + dy[i], now.dis + 1, i};
        if (check_range(nxt.x, nxt.y) && check_map(nxt.x, nxt.y) &&
            !vis[nxt.x][nxt.y]) {
          q.push(nxt);
          to_berth[nxt.x][nxt.y].push_back({berth_id, i ^ 1, nxt.dis});
          vis[nxt.x][nxt.y] = 1;
        }
      }
    }
  }
} berth[berth_num + 10];

struct Boat {
  int num, pos, status;
} boat[10];

void Init() {
  freopen("test.in", "r", stdin);
  freopen("test.out", "w", stdout);
  for (int i = 0; i < n; i++) scanf("%s", ch[i]);
  for (int i = 0; i < berth_num; i++) {
    int id;
    scanf("%d", &id);
    berth[id].berth_id = id;
    scanf("%d%d%d%d", &berth[id].x, &berth[id].y, &berth[id].transport_time,
          &berth[id].loading_speed);
  }
  for (int i = 0; i < berth_num; ++i) {
    berth[i].berth_bfs();
  }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (check_map(i, j)) sort(to_berth[i][j].begin(), to_berth[i][j].end());
  scanf("%d", &boat_capacity);
  char okk[100];
  scanf("%s", okk);
  printf("OK\n");
  fflush(stdout);
}

class mcmf {
 private:
  struct Edge {
    int vet, nxt, f;
    double val;
  } e[N * N];
  int head[N * N], tot;

 public:
  void add(int u, int v, int f, double val) {
    e[++tot] = {v, head[u], f, val}, head[u] = tot;
    e[++tot] = {u, head[v], 0, -val}, head[v] = tot;
  }
};

int Input() {
  scanf("%d%d", &frame_id, &money);
  int num;
  scanf("%d", &num);
  for (int i = 1; i <= num; i++) {
    int x, y, val;
    scanf("%d%d%d", &x, &y, &val);
  }
  for (int i = 0; i < robot_num; i++) {
    int sts;
    scanf("%d%d%d%d", &robot[i].goods, &robot[i].x, &robot[i].y, &sts);
  }
  for (int i = 0; i < 5; i++) scanf("%d%d\n", &boat[i].status, &boat[i].pos);
  char okk[100];
  scanf("%s", okk);
  return frame_id;
}

int main() {
  Init();
  return 0;
  for (int zhen = 1; zhen <= 15000; zhen++) {
    int id = Input();
    mcmf();
    for (int i = 0; i < robot_num; i++) printf("move %d %d\n", i, rand() % 4);
    puts("OK");
    fflush(stdout);
  }

  return 0;
}
