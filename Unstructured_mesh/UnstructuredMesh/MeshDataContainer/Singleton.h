#ifndef SINGLETON_H
#define SINGLETON_H
#include <type_traits>
#include <memory>


/**
 * @brief The Singleton class
 */
template<class Class>
class Singleton {
    //static_assert (std::is_trivially_constructible<Class>::value, "The class in singleton must be trivially constructible.");
public:
    static Class& getInstance(){
        if (p == nullptr) {
            p = std::unique_ptr<Class>(new Class); // the class must be trivially constructible
        }
        return *p;
    }
protected:
    static std::unique_ptr<Class> p;
    Singleton() {}
private:
    // disable move and copy options
    Singleton(Singleton const&) = delete;
    Singleton(Singleton const&&) = delete;
    Singleton& operator=(Singleton const&) = delete;
    Singleton& operator=(Singleton const&&) = delete;

};
template <class Class> std::unique_ptr<Class> Singleton<Class>::p = nullptr;



#endif // SINGLETON_H
