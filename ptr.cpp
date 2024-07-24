#include <iostream>

template <typename T>
class sharedptr{
public:
    sharedptr(T* ptr = nullptr) : ptr(ptr), count(new int(1)){}

    sharedptr(const sharedptr& pother) : ptr(pother.ptr), count(pother.count){
        ++*count;
    }

    sharedptr& operator=(const sharedptr& pother){
        if(ptr == pother.ptr){
            return *this;
        }
        if(ptr){
            release();
            ptr = pother.ptr;
            count = pother.count;
            ++*count;
        }
        return *this;
    }

    ~sharedptr(){
        release();
    }

    T& operator*(){
        return * ptr;
    }

    T* operator->(){
        return ptr;
    }
    int use_count(){
        return *count;
    }
private:
    void release(){
        if(--*count == 0){
            delete ptr;
            delete count;
        }
    }
    T* ptr;
    int* count;
};

class MyClass {
public:
    MyClass() { std::cout << "MyClass 构造函数\n"; }
    ~MyClass() { std::cout << "MyClass 析构函数\n"; }
    void do_something() { std::cout << "MyClass::do_something() 被调用\n"; }
};


int main() {
    {
        sharedptr<MyClass> ptr1(new MyClass());
        {
            sharedptr<MyClass> ptr2 = ptr1;
            ptr1->do_something();
            ptr2->do_something();
            std::cout << "引用计数: " << ptr1.use_count() << std::endl;
        }
        std::cout << "引用计数: " << ptr1.use_count() << std::endl;
    }

    return 0;
}
