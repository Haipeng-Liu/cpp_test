#include <bits/stdc++.h>
//#include "kalman.h"
#include "math_test.h"
#include <algorithm>
#include <thread>
#include "thread_test.h"

//using namespace std;


struct PoseInfo {
    float x;
    float y;
};

typedef PoseInfo Point2D;
typedef PoseInfo Vector2D;
float Vector2DAngle(const Vector2D& vec1, const Vector2D& vec2)
{
//    double PI = 3.141592653;
    float t = (vec1.x * vec2.x + vec1.y * vec2.y) / (sqrt(pow(vec1.x, 2) + pow(vec1.y, 2)) * sqrt(pow(vec2.x, 2) + pow(vec2.y, 2)));
    float angle = acos(t) * (180 / PI);
    return angle;
}

void testVector2DAngle()
{
    Point2D vecA_start_point;
    vecA_start_point.x = 0.0;
    vecA_start_point.y = 0.0;

    Point2D vecA_end_point;
    vecA_end_point.x = 1.0;
    vecA_end_point.y = 0.0;

    Point2D vecB_end_point;
    vecB_end_point.x = 1.0;
    vecB_end_point.y = 5.0;

    Vector2D vecA;
    vecA.x = vecA_end_point.x - vecA_start_point.x;
    vecA.y = vecA_end_point.y - vecA_start_point.y;

    Vector2D vecB;
    vecB.x = vecB_end_point.x - vecA_start_point.x;
    vecB.y = vecB_end_point.y - vecA_start_point.y;

    float angle = Vector2DAngle(vecA, vecB);

    std::cout << "angle = " << angle << std::endl;
}

// 定义向量结构体
struct Vector {
    double x;
    double y;
    double z;
};

// 计算两个向量的夹角（弧度制）
double angle(Vector v1, Vector v2) {
    // 计算两个向量的点积
    double dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

    // 计算两个向量的模长
    double len1 = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
    double len2 = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);

    // 计算两个向量的夹角（弧度制）
    double angle = acos(dot / (len1 * len2));

    return angle;
}

void angleTest()
{
    // 定义两个向量
    Vector v1 = {1, 2, 3};
    Vector v2 = {4, 5, 6};

    // 计算两个向量的夹角（弧度制）
    double angleRad = angle(v1, v2);

    // 将弧度转换为角度
    double angleDeg = angleRad * 180 / M_PI;

    std::cout << "The angle between the two vectors is " << angleDeg << " degrees." << std::endl;
}

const int MAX_ITERATIONS = 10000;
const double ACCURACY = 0.00001;

void get_eigenvectors(double matrix[3][3], double *eigenvalues, double *eigenvectors)
{
    double b[3] = {1, 1, 1};
    double x[3];
    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        double s = 0;
        for (int j = 0; j < 3; j++)
        {
            x[j] = 0;
            for (int k = 0; k < 3; k++)
            {
                x[j] += matrix[j][k] * b[k];
            }
            s += x[j] * x[j];
        }

        s = sqrt(s);
        for (int j = 0; j < 3; j++)
        {
            b[j] = x[j] / s;
        }

        double error = 0;
        for (int j = 0; j < 3; j++)
        {
            error += std::abs(x[j] - b[j]);
        }

        if (error < ACCURACY)
        {
            break;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        eigenvalues[i] = 0;
        for (int j = 0; j < 3; j++)
        {
            eigenvalues[i] += matrix[i][j] * b[j];
        }
        eigenvectors[i] = b[i];
    }
}

int getUssObjPos(double x1, double y1, double r1,
                 double x2, double y2, double r2)
{
//    double x1, y1, r1, x2, y2, r2;
//    // 输入第一个圆的数据（x坐标、y坐标和半径）
////    cin >> x1 >> y1 >> r1;
//    x1 = 1; y1 = 0; r1 = 1;
//    // 输入第二个圆的数据（x坐标、y坐标和半径）
////    cin >> x2 >> y2 >> r2;
//    x2 = 3; y2 = 0; r2 = 2;

    // 计算两个圆心之间的距离
    double d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

    // 如果两圆相离或包含，则无交点
    if (d > r1 + r2 || d < std::abs(r1 - r2)) {
        std::cout << "No intersection point." << std::endl;
        return 0;
    }

    // 计算两个圆的交点
    double a = (pow(r1, 2) - pow(r2, 2) + pow(d, 2)) / (2 * d);
    double h = sqrt(pow(r1, 2) - pow(a, 2));
    double cx = x1 + a * (x2 - x1) / d;
    double cy = y1 + a * (y2 - y1) / d;
    double ix1 = cx + h * (y2 - y1) / d;
    double ix2 = cx - h * (y2 - y1) / d;
    double iy1 = cy - h * (x2 - x1) / d;
    double iy2 = cy + h * (x2 - x1) / d;

    // 输出交点的坐标
    std::cout << "(" << ix1 << ", " << iy1 << ")" << std::endl;
    std::cout << "(" << ix2 << ", " << iy2 << ")" << std::endl;

    return 0;
}

void findIntersection(double x1, double y1, double r1,
                      double x2, double y2, double r2) {
    // 计算���圆心之间的距离
    double d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

    // 如果两圆相离或包含，则无交点
    if (d > r1 + r2 || d < abs(r1 - r2)) {
        std::cout << "No intersection point." << std::endl;
        return;
    }

    // 计算交点的坐标
    double a = (pow(r1, 2) - pow(r2, 2) + pow(d, 2)) / (2 * d);
    double h = sqrt(pow(r1, 2) - pow(a, 2));
    double x3 = x1 + a * (x2 - x1) / d;
    double y3 = y1 + a * (y2 - y1) / d;
    double x4 = x3 + h * (y2 - y1) / d;
    double y4 = y3 - h * (x2 - x1) / d;
    double x5 = x3 - h * (y2 - y1) / d;
    double y5 = y3 + h * (x2 - x1) / d;

    // 输出交点坐标
    std::cout << "Intersection points: (" << x4 << ", " << y4 << ") and ("
         << x5 << ", " << y5 << ")" << std::endl;

    // 计算交点的坐标
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
    h = sqrt(r1 * r1 - a * a);
    x3 = x1 + a * (x2 - x1) / d;
    y3 = y1 + a * (y2 - y1) / d;
    x4 = x3 + h * (y2 - y1) / d;
    y4 = y3 - h * (x2 - x1) / d;

    // 输出交点坐标
    std::cout << "(" << x4 << "," << y4 << ")" << std::endl;

    float a1, b1, R1, a2, b2, R2;
    a1 = x1;
    b1 = y1;
    R1 = r1;

    a2 = x2;
    b2 = y2;
    R2 = r2;

    //
    float R1R1 = R1*R1;
    float a1a1 = a1*a1;
    float b1b1 = b1*b1;

    float a2a2 = a2*a2;
    float b2b2 = b2*b2;
    float R2R2 = R2*R2;

    float subs1 = a1a1 - 2 * a1*a2 + a2a2 + b1b1 - 2 * b1*b2 + b2b2;
    float subs2 = -R1R1 * a1 + R1R1 * a2 + R2R2 * a1 - R2R2 * a2 + a1a1*a1 - a1a1 * a2 - a1*a2a2 + a1*b1b1 - 2 * a1*b1*b2 + a1*b2b2 + a2a2*a2 + a2*b1b1 - 2 * a2*b1*b2 + a2*b2b2;
    float subs3 = -R1R1 * b1 + R1R1 * b2 + R2R2 * b1 - R2R2 * b2 + a1a1*b1 + a1a1 * b2 - 2 * a1*a2*b1 - 2 * a1*a2*b2 + a2a2 * b1 + a2a2 * b2 + b1b1*b1 - b1b1 * b2 - b1*b2b2 + b2b2*b2;
    float sigma = sqrt((R1R1 + 2 * R1*R2 + R2R2 - a1a1 + 2 * a1*a2 - a2a2 - b1b1 + 2 * b1*b2 - b2b2)*(-R1R1 + 2 * R1*R2 - R2R2 + subs1));
    if(abs(subs1)>0.0000001)//分母不为0
    {
        x4 = (subs2 - sigma*b1 + sigma*b2) / (2 * subs1);
        x5 = (subs2 + sigma*b1 - sigma*b2) / (2 * subs1);

        y4 = (subs3 + sigma*a1 - sigma*a2) / (2 * subs1);
        y5 = (subs3 - sigma*a1 + sigma*a2) / (2 * subs1);
    }
    std::cout << "Intersection points: (" << x4 << ", " << y4 << ") and ("
         << x5 << ", " << y5 << ")" << std::endl;
}

void testFindIntersection()
{
    double x1, y1, r1, x2, y2, r2;
    // 输入第一个圆的数据（x坐标、y坐标和半径）
    x1 = 1; y1 = 0; r1 = 1;
    // 输入第二个圆的数据（x坐标、y坐标和半径）
    x2 = 3; y2 = 0; r2 = 2;
    std::cout << "test 1: " << std::endl;
//    getUssObjPos(x1, y1, r1, x2, y2, r2);
    findIntersection(x1, y1, r1, x2, y2, r2);

    x1 = 0; y1 = 0; r1 = 1;
    x2 = 5; y2 = 5; r2 = 2;
    std::cout << "test 2: " << std::endl;
//    getUssObjPos(x1, y1, r1, x2, y2, r2);
    findIntersection(x1, y1, r1, x2, y2, r2);

    std::cout << "test 3: " << std::endl;
    x1 = 0; y1 = 0; r1 = 1;
    x2 = 3; y2 = 0; r2 = 2;
//    getUssObjPos(x1, y1, r1, x2, y2, r2);
    findIntersection(x1, y1, r1, x2, y2, r2);

    std::cout << "test 4: " << std::endl;
    x1 = 0; y1 = 0; r1 = 1;
    x2 = 1; y2 = 0; r2 = 2;
//    getUssObjPos(x1, y1, r1, x2, y2, r2);
    findIntersection(x1, y1, r1, x2, y2, r2);

    std::cout << "test 5: " << std::endl;
    x1 = 0; y1 = 0; r1 = 2;
    x2 = 0; y2 = 0; r2 = 1;
//    getUssObjPos(x1, y1, r1, x2, y2, r2);
    findIntersection(x1, y1, r1, x2, y2, r2);
}

void test_get_eigenvectors()
{

//    double matrix[3][3] = {{1, 2, 3},
//                           {2, 4, 5},
//                           {3, 5, 6}};
    double matrix[3][3] = {{3, 2, 4},
                           {2, 0, 2},
                           {4, 2, 3}};

    double eigenvalues[3];
    double eigenvectors[3];
    get_eigenvectors(matrix, eigenvalues, eigenvectors);

    for (int i = 0; i < 3; i++)
    {
        std::cout << "Eigenvalue: " << eigenvalues[i] << std::endl;
        std::cout << "Eigenvector: (" << eigenvectors[i] << ", " << eigenvectors[i + 1] << ", " << eigenvectors[i + 2] << ")" << std::endl;
    }
}

struct studentInfo {
    int index;
    std::string name;
    int price;
    int develop;
};

void questionC()
{
    int num = 0;
    std::cin >> num;
    std::vector<studentInfo> students;
    for (int i = 0; i < num; i++) {
        studentInfo info;
        info.index = i+1;
//        std::cin >> info.name;
//        getchar();
//        std::cin >> info.price;
//        getchar();
//        std::cin >> info.develop;
//        getchar();

        std::cin >> info.name >> info.price >> info.develop;

        students.emplace_back(info);
    }

    std::sort(students.begin(), students.end(), [](studentInfo& a, studentInfo& b) {
       if (a.price != b.price) return a.price > b.price;
       else if (a.develop != b.develop) return a.develop > b.develop;
       else return a.index < b.index;
    });
    for (auto& item : students) {
        std::cout << item.index << " " << item.name << " " << item.price << " " << item.develop << " " << item.price * item.develop << std::endl;
    }
}

int questionD()
{
    int count = 0;
    std::cin >> count;
    if (count == 0) {
        return 0;
    }
    std::priority_queue<int, std::vector<int>, std::greater<int>> que;
    for (int i = 0; i < count; i++) {
        int weight = 0;
        std::cin >> weight;
        que.emplace(weight);
    }
    int ans = 0;
    while(!que.empty()) {
        int a = que.top();
        que.pop();
        int b = que.top();
        que.pop();
        ans += (a + b);
        que.emplace(a + b);
    }
    std::cout << ans << std::endl;
    return ans;
}

int questionE()
{
    int count = 0;
    std::cin >> count;
    if (count == 0) {
        return 0;
    }
    std::priority_queue<long, std::vector<long>, std::greater<long>> que;
    for (int i = 0; i < count; i++) {
        int weight = 0;
        std::cin >> weight;
        que.emplace(weight);
    }
    long ans = 0;
    while(que.size() >= 3) {
        long a = que.top();
        que.pop();
        long b = que.top();
        que.pop();
        long c = que.top();
        que.pop();
        ans += (a + b + c);
        que.emplace(a + b + c);
    }
    if (que.size() == 1)
        std::cout << ans << std::endl;
    else {
        std::cout << "9223372036854775807" << std::endl;
    }
    return 0;
}

const int MAXN = 100005;

int w[MAXN];
int head[MAXN], nxt[MAXN << 1], to[MAXN << 1], cnt;

void addEdge(int u, int v) {
    nxt[++cnt] = head[u];
    head[u] = cnt;
    to[cnt] = v;
}

int f[MAXN][2];

void dfs(int u, int fa) {
    f[u][0] = 0;
    f[u][1] = w[u];
    for (int i = head[u]; i; i = nxt[i]) {
        int v = to[i];
        if (v == fa) continue;
        dfs(v, u);
        f[u][0] += std::max(f[v][0], f[v][1]);
        f[u][1] += f[v][0];
    }
}

int questionF()
{
    int n;
    scanf("%d", &n);
    for (int i = 1; i <= n; i++) {
        scanf("%d", &w[i]);
    }
    for (int i = 1; i < n; i++) {
        int u, v;
        scanf("%d%d", &u, &v);
        addEdge(u, v);
        addEdge(v, u);
    }
    dfs(1, 0);
    printf("%d\n", std::max(f[1][0], f[1][1]));
    return 0;

    int count = 0;
    std::cin >> count;
    if (count == 0) {
        return 0;
    }
    int ans = 0;
    std::map<int, std::set<int>> nodeMap;
    for (int i = 0; i < count; i++) {
        int weight = 0;
        std::cin >> weight;
    }
    for (int i = 0; i < count-1; i++) {
        int a,b;
        std::cin >> a >> b;
        nodeMap[a].emplace(b);
        nodeMap[b].emplace(a);
    }
    for (int i = count; i > 0 ; i--) {
        if (!nodeMap.count(i)) {
            continue;
        }
        ans += i;
        for (auto a : nodeMap[i]) {
            nodeMap.erase(a);
        }
    }

    std::cout << ans << std::endl;
    return 0;
}

// 使用变参模板实现递归打印函数
void print() {
    std::cout << std::endl;
}
template<typename T, typename... Args>
void print(T first, Args... args) {
    std::cout << first << " ";
    print(args...);
}

int main()
{
//    testVector2DAngle();


    std::string str2{"123\0abc", 5};
    uint8_t cc = 56; // 这里是单纯想放一个整数
    std::cout << cc << std::endl; // 但这里会打印出8，而不是56

    std::vector<int> tempVt;
    std::cout << "sizeof(std::vector<int>) :" << sizeof(std::vector<int>) << " sizeof(int):" << sizeof(int) << ", size " << tempVt.size() << ", " << tempVt.capacity() << std::endl;   // output: 24 4


//    test_distSegmentToRectangle();
//    test_Angle();
//
//    questionF();



//    print(1, 2, 3, "Hello", 4.5);  // 调用变参模板函数 print


    ThreadTest::threadTest1();

    return 0;
}

std::mutex m;
std::condition_variable cv;
std::string data_str;
bool ready = false;
bool processed = false;

void worker_thread()
{
    std::cout << "3、worker_thread子线程开始执行"  << std::endl;
    // Wait until main() sends data
    std::unique_lock<std::mutex> lk(m);
    std::cout << "4、worker_thread子线程获取到锁,条件满足无需notify，不阻塞向下执行"  << std::endl;
    cv.wait(lk, []{return ready;});

    // after the wait, we own the lock.
    data_str += " after processing";
    // Send data back to main()
    processed = true;
    std::cout << "5、Worker thread signals data processing completed\n";

    // Manual unlocking is done before notifying, to avoid waking up
    // the waiting thread only to block again (see notify_one for details)
    lk.unlock();
    std::cout << "6、worker_thread子线程交出执行权限，主线程执行"  << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    cv.notify_one();
    std::cout << "9、worker_thread调用 notify_one"  << std::endl;
}
int threadTest()
{
    std::thread worker(worker_thread);
    std::cout << "1、主线程开始执行"  << std::endl;
    data_str = "Example data";
    // send data to the worker thread
    {
        //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        std::lock_guard<std::mutex> lk(m);
        ready = true;
    }
    std::cout << "2、锁已经释放了，主线程休眠，子线程执行"  << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    //cv.notify_one();
    {
        std::cout << "7、主线程data：" << data_str << std::endl;
        std::unique_lock<std::mutex> lk(m);
        std::cout << "8、主线程条件满足无需notify" << std::endl;
        cv.wait(lk, []{return processed;});
    }

    worker.join();
    std::cout << "10、主线程结束" << std::endl;

    return 0;
}
