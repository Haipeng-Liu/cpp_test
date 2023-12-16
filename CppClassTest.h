//
// Created by 刘海鹏 on 2023/12/16.
//

#ifndef TEST_CPPCLASSTEST_H
#define TEST_CPPCLASSTEST_H

#include <map>
#include <unordered_map>

using namespace std;
class CppClassTest {
public:
    static int add_func(int a, int b)
    {
        return a + b;
    }
    int multiply_func(int a, int b)
    {
        return a * b;
    }

    void testFunction()
    {
        typedef int(*Func)(int ,int);
        Func f1 = add_func;
        cout<< f1(1,2)<<endl;  //3

        std::function<int(int, int)> f2 = add_func;
        cout<<f2(1, 2)<<endl;      // 3
    }

    void Function1() {
        func_value_++;
        std::cout << "Executing Function 1, value " << func_value_ << std::endl;
    }

    void Function2() {
        func_value_++;
        std::cout << "Executing Function 2, value " << func_value_ << std::endl;
    }

    void Function3() {
        func_value_++;
        std::cout << "Executing Function 3, value " << func_value_ << std::endl;
    }
    // 定义一个函数指针类型，指向类内的成员函数
    typedef void (CppClassTest::*FunctionPointer)();
    // 使用std::map存储key和函数指针的映射关系
    std::map<int, FunctionPointer> functionMap = {
            {1, &CppClassTest::Function1},
            {2, &CppClassTest::Function2},
            {3, &CppClassTest::Function3}
    };

    // 根据key值执行相应的函数
    void ExecuteFunction(int key) {
        if (functionMap.find(key) != functionMap.end()) {
            (this->*functionMap[key])();
        } else {
            std::cout << "Key not found" << std::endl;
        }
    }

private:
    int func_value_ = 0;
    std::unordered_map<uint8_t, std::function<int(int, int)>> func_map;
};
#endif //TEST_CPPCLASSTEST_H
