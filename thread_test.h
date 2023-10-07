//
// Created by dell on 2023/8/29.
//

#ifndef TEST_THREAD_TEST_H
#define TEST_THREAD_TEST_H
#include <iostream>
#include <thread>       // std::thread, std::this_thread::yield
#include <mutex>          // std::mutex, std::unique_lock
#include <atomic>
#include <condition_variable> // std::condition_variable

using namespace std;

namespace ThreadTest {

void thread_1()
{
    while(1)
    {
        cout<<"子线程1111"<<endl;
    }
}

void thread_2(int x)
{
    while(1)
    {
        cout<<"子线程2222"<<endl;
    }
}

std::mutex mtx;           // mutex for critical section
void print_block (int n, char c) {
    // critical section (exclusive access to std::cout signaled by locking mtx):
    mtx.lock();
    for (int i=0; i<n; ++i) { std::cout << c; }
    std::cout << '\n';
    mtx.unlock();
}

int g_i = 0;
std::mutex g_i_mutex;  // protects g_i，用来保护g_i
void safe_increment()
{
    const std::lock_guard<std::mutex> lock(g_i_mutex);
    ++g_i;
    std::cout << std::this_thread::get_id() << ": " << g_i << '\n';
    // g_i_mutex自动解锁
}


struct Box {
    explicit Box(int num) : num_things{num} {}

    int num_things;
    std::mutex m;
};
void transfer(Box &from, Box &to, int num)
{
    // defer_lock表示暂时unlock，默认自动加锁
    std::unique_lock<std::mutex> lock1(from.m, std::defer_lock);
    std::unique_lock<std::mutex> lock2(to.m, std::defer_lock);

    //两个同时加锁
    std::lock(lock1, lock2);//或者使用lock1.lock()

    from.num_things -= num;
    to.num_things += num;
    //作用域结束自动解锁,也可以使用lock1.unlock()手动解锁
}

atomic_int n1(0);//std::atomic_int只是std::atomic<int>的别名罢了。
void count10000() {
    for (int i = 1; i <= 10000; i++) {
        n1++;
    }
}


std::condition_variable cv;
int cargo = 0;
bool shipment_available() {return cargo!=0;}
void consume (int n) {
    for (int i=0; i<n; ++i) {
        std::unique_lock<std::mutex> lck(mtx);//自动上锁
        //第二个参数为false才阻塞（wait），阻塞完即unlock，给其它线程资源
        cv.wait(lck,shipment_available);
        // consume:
        std::cout << cargo << '\n';
        cargo=0;
    }
}


int value;
void read_value() {
    std::cin >> value;
    cv.notify_one();
}


int mythread() //线程入口函数
{
    cout << "mythread start" << " threadid= " << std::this_thread::get_id() << endl; //打印线程id

    std::chrono::milliseconds dura(5000); //定一个5秒的时间
    std::this_thread::sleep_for(dura);  //休息一定时常

    cout << "mythread end" << " threadid= " << std::this_thread::get_id() << endl; //打印线程id

    return 5;
}

template<class ... Args> decltype(auto) sum(Args&&... args) {
    return (0 + ... + args);
}

template<class ... Args> void sum_thread(promise<long long> &val, Args&&... args) {
    val.set_value(sum(args...));
}


int threadTest1()
{

//    thread first ( thread_1);  // 开启线程，调用：thread_1()
//    thread second (thread_2,100); // 开启线程，调用：thread_2(100)
//    first.detach();
//    second.detach();
//    for(int i = 0; i < 10; i++)
//    {
//        std::cout << "主线程\n";
//    }

//    std::thread th1 (print_block,50,'*');//线程1：打印*
//    std::thread th2 (print_block,50,'$');//线程2：打印$
//    th1.join();
//    th2.join();

//    std::cout << "main id: " <<std::this_thread::get_id()<<std::endl;
//    std::cout << "main: " << g_i << '\n';
//    std::thread t1(safe_increment);
//    std::thread t2(safe_increment);
//    t1.join();
//    t2.join();
//    std::cout << "main: " << g_i << '\n';

//    Box acc1(100);
//    Box acc2(50);
//    std::thread t1(transfer, std::ref(acc1), std::ref(acc2), 10);
//    std::thread t2(transfer, std::ref(acc2), std::ref(acc1), 5);
//    t1.join();
//    t2.join();
//    std::cout << "acc1 num_things: " << acc1.num_things << std::endl;
//    std::cout << "acc2 num_things: " << acc2.num_things << std::endl;

//    thread th[100];
//    for (thread &x : th)
//        x = thread(count10000);
//    for (thread &x : th)
//        x.join();
//    cout << n << endl;

//    std::thread consumer_thread (consume,10);
//    for (int i=0; i<10; ++i) {
//        //每次cargo每次为0才运行。
//        while (shipment_available()) std::this_thread::yield();
//        std::unique_lock<std::mutex> lck(mtx);
//        cargo = i+1;
//        cv.notify_one();
//    }
//    consumer_thread.join();


//    std::cout << "Please, enter an integer (I'll be printing dots): \n";
//    std::thread th (read_value);
//    std::mutex mtx;
//    std::unique_lock<std::mutex> lck(mtx);
//    while (cv.wait_for(lck,std::chrono::seconds(1))==std::cv_status::timeout) {
//        std::cout << '.' << std::endl;
//    }
//    std::cout << "You entered: " << value << '\n';
//    th.join();


    /************************ std::async ***************/
//    cout << "main " << "threadid= " << std::this_thread::get_id() << endl;
//    std::future<int> result = std::async(mythread);//流程并不卡在这里
//    cout << "continue....." << endl;
//    //枚举类型
//    std::future_status status = result.wait_for(std::chrono::seconds(0));//等待一秒
//    if (status == std::future_status::deferred)
//    {
//        //线程被延迟执行了，系统资源紧张
//        cout << result.get() << endl; //此时采取调用mythread()
//    }
//    else if (status == std::future_status::timeout)//
//    {
//        //超时：表示线程还没执行完；我想等待你1秒，希望你返回，你没有返回，那么 status = timeout
//        //线程还没执行完
//        cout << "超时：表示线程还没执行完!" << endl;
//    }
//    else if (status == std::future_status::ready)
//    {
//        //表示线程成功返回
//        cout << "线程成功执行完毕，返回!" << endl;
//        cout << result.get() << endl;
//    }
//    cout << "I love China!" << endl;

//    async(launch::async, [](const char *message){
//        cout << message << flush;
//    }, "Hello, ");
//    cout << "World!" << endl;


    /************************ std::future ***************/
    //通过async来获取异步操作结果
//    std::future<int>  result = std::async([](){
//        std::this_thread::sleep_for(std::chrono::milliseconds(500));
//        return 8;
//    });
//    std::cout << "the future result : " << result.get() << std::endl;
//    std::cout << "the future status : " << result.valid() << std::endl;
//    try
//    {
//        result.wait();  //或者 result.get() ,会异常
//        //因此std::future只能用于单线程中调用 ，多线程调用使用std::shared_future();
//    }
//    catch (...)
//    {
//        std::cout << "get error....\n ";
//    }

/************************ std::promise ***************/
    promise<long long> sum_value;
    thread get_sum(sum_thread<int, int, int>, ref(sum_value), 1, 10, 100);
    cout << sum_value.get_future().get() << endl;
    get_sum.join();


/************************ std::packaged_task ***************/
    std::packaged_task<int()> task([](){return 7;});
    std::thread t1(std::ref(task));
    std::future<int> f1 = task.get_future();
    auto r1 = f1.get();
    cout << r1 << endl;


    return 0;
}

}

#endif //TEST_THREAD_TEST_H
