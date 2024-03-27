#include <bits/stdc++.h>
#define PFI pair<float, int>
const float k1, k2, k3, k4;
const int dividing_frame = 14995;
int frame_time, money;
#include "MyMap.hpp"
#include "Robot_Boat_Distribute.hpp"

const int n = 200, robot_num = 10, berth_num = 10, boat_num = 5;

void Init() { MyMap &mp = MyMap::GetMap(); }

using namespace std;
int main() {
  Init();
  Robot_Boat_Distribute RBD(n, robot_num, berth_num, boat_num);
  RBD.Init();
  for (int zhen = 1; zhen <= 15000; ++zhen) {
    scanf("%d%d", &frame_time, &money);
    RBD.frame_input();
  }
  delete (&MyMap::GetMap());
  delete (&Goods_Table::GetGT());
  return 0;
}