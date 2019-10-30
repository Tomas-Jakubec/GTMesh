#include "../debug/Debug.h"
#include "../Unstructured_mesh/UnstructuredMesh/UnstructuredMesh.h"
#include "../Unstructured_mesh/UnstructuredMesh/MeshDataContainer/MemberApproach.h"
#include "../Unstructured_mesh/UnstructuredMesh/MeshDataContainer/Traits.h"
#include <functional>
#include <type_traits>
#include <iostream>
#include <list>
#include <map>
using namespace std;


template <typename T1, typename T2 = void>
struct has_public_member : public std::integral_constant<bool, false> {

};

template <typename T1>
struct has_public_member<T1, typename std::enable_if<
        noexcept (T1::data) || !noexcept (T1::data)>::type> : public std::integral_constant<bool, true> {

};

template <typename ResulType, typename ...ArgTypes>
struct function_ptr {
    typedef ResulType(*type)(ArgTypes...);
};


template <typename Class,typename ResulType, typename ...ArgTypes>
struct member_function_ptr {
    using type = ResulType (Class::*) (ArgTypes...);
};

template <typename Class,typename ResulType, typename ...ArgTypes>
struct member_const_function_ptr {
    using type = ResulType (Class::*) (ArgTypes...) const;
};


template <typename Class,typename MemberType>
struct member_ptr {
    typedef MemberType Class::* type;
};



struct Temp {
    double data;
    double fun(double d){return data;}
};

void testDebug() {
    double r = 42.15;
    int i = 15;
    char c = 42;
    bool b = false;
    std::list<int> list = {1,2,3};
    std::vector<std::list<int>> vec(5, list);
    std::map<std::string, size_t> m{
        {"prvni", 1},
        {"druhy", 2},
        {"treti", 3}
    };
    ConsoleLogger::writeVar(__LINE__, __FILE__, "r", r, "i", i, "c", c, "list", list, "vec", vec, "b", b, "map", m);
    ConsoleLogger::writeVar(__LINE__, __FILE__,"---", {5,4,3,2});
    DBGVAR(r, i, c, list, vec, b, m);

    Vertex<7, double> vert;
    DBGVAR(vert, vert.getCoordinates());

    DBGVAR((Detail::is_exportable<decltype(vec)>::value));

    DBGVAR((Detail::is_exportable<double>::value));

    DBGVAR(Detail::is_indexable<double>::value);

    DBGVAR(Detail::is_indexable<decltype(vec)>::value);

    DBGVAR(Detail::is_indexable<decltype(vert)>::value);

    Subelement<size_t> s({1,true});
    DBGVAR(s);

    HTMLDBGVAR(r, i, c, list, vec, b, m);

    HTMLDBGVAR(r+1, i+1, char(c+1), list, vec[0], b, m["prvni"]);
}



//The concept implementation
template<typename T>
class NeedIterator{
    static_assert (Detail::is_iterable<T>::value, "The type must be iterable");
public:
    NeedIterator(const T&){}
};

template <unsigned int ... Is>
class ClassA
 {
   public:
      ClassA (std::integer_sequence<unsigned int,Is...>)
      {DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...})); }

      static void fun (std::index_sequence<Is...>)
      {DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...})); }

 };


//template <typename ... t> class ClassB{};

template<typename  Tuple, unsigned int ... Is>
class ClassB
 {
   public:
      ClassB (const std::integer_sequence<unsigned int, Is...>, Tuple t)
      {std::tuple_element_t<0,Tuple> typ = 0;
          DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...}), std::get<0>(t), typ); }

 };


template <typename ... t> class ClassC{};

template<unsigned int ... Is, typename ...Types>
class ClassC<std::integer_sequence<unsigned int,Is...>, std::tuple<Types...>>
 {
   public:
      ClassC (const std::integer_sequence<unsigned int, Is...>, std::tuple<Types...>)
      {std::tuple_element_t<0,std::tuple<Types...>> typ = 0;
          DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...}), typ); }

      ClassC () {
          std::tuple_element_t<0,std::tuple<Types...>> typ = 42.15;
          std::tuple_element_t<1,std::tuple<Types...>> typ2 = 42.15;
          DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...}), typ, typ2);
      }
 };



void testTemplate() {
    ClassA n(std::make_integer_sequence<unsigned int, 3>{});
    UnstructuredMesh<3,size_t, double,6> mesh3;
    //MeshDataContainer<Vertex<3, double>, 0,1,2> centers2(mesh3,std::make_integer_sequence<unsigned int, 3>{}, Vertex<3, double>{});
    //ComputeCenters(mesh3);

    ClassA p(make_custom_integer_sequence_t<unsigned int, 10, 0, -2>{});
    std::tuple<double, char> t{};
    t={1,2};
    ClassB u(make_custom_integer_sequence_t<unsigned int, 2, 0, -2>{}, t);

    ClassC<std::integer_sequence<unsigned int, 2,0>, std::tuple<double, char>> c(make_custom_integer_sequence_t<unsigned int, 2, 0, -2>{}, std::tuple<double, char>{});
    ClassC<std::integer_sequence<unsigned int, 2,0>, std::tuple<double, char>> cc;
    ClassC<std::integer_sequence<unsigned int, 2,0>, decltype(std::make_tuple(1.0, 'a'))> ccc;

    NeedIterator valid(mesh3.getCells());
    //NeedIterator invalid(0.0);
}


template<typename T, typename VOID = void>
struct is_indexable : public integral_constant<bool, false>{};

template<typename T>
struct is_indexable<
        T,
        enable_if<
        is_same<decltype(&T::operator[]), typename member_const_function_ptr<vector<double>,const double&, size_t>::type>::value
        >

        > : public integral_constant<bool, false>{};


template <typename T>
class OutOfLineSpecialization {
public:
    template<typename ... Sequence>
    static void printSeq(Sequence...);

    template<T... data>
    static void onlySpec();
};

template<typename T> template<typename...Sequence>
void OutOfLineSpecialization<T>::printSeq(Sequence... onlyOne) {
    DBGMSG("capture");
    DBGVAR(onlyOne...);
}
template<> template<typename...Sequence>
void OutOfLineSpecialization<int>::printSeq(Sequence... onlyOne) {
    DBGMSG("int");
    DBGVAR(onlyOne...);
}


template<> template<int...data>
void OutOfLineSpecialization<int>::onlySpec() {
    DBGMSG("int");
    DBGVAR(data...);
}



template <typename T>
class InLineSpecialization {
public:
    template<typename...Sequence>
    static void printSeq(Sequence... onlyOne) {
        DBGMSG("capture");
        DBGVAR(onlyOne...);
    }
};


template <>
class InLineSpecialization<int> {
public:
    template<typename...Sequence>
    static void printSeq(Sequence... onlyOne) {
        DBGMSG("int");
        DBGVAR(onlyOne...);
    }
    template<int...data>
    static void onlySpec() {
        DBGMSG("int");
        DBGVAR(data...);
    }
};


template<typename T>
void callSpec() {
    // supr, tohle nabízí
    OutOfLineSpecialization<T>::printSeq(1);
    OutOfLineSpecialization<T>::template onlySpec<2>();

    InLineSpecialization<T>::printSeq(1);
    InLineSpecialization<T>::template onlySpec<3>();
}


struct tempData {
    double density;

    Vector<3,double> velocity;

    double& getData(){
        return density;
    }

    Vector<3,double> getMomentum(){
        return velocity*density;
    }

    void setMomentum(const Vector<3,double>& val){
        velocity = val / density;
    }

};


void testMemberRef(){

    MemberApproach<tempData, double>* app;

    //MemberReference<tempData, double, bool> invalid;

    MemberReference<tempData, double>::SuperRef<double tempData::*> ref(&tempData::density);
    MemberReference<tempData, double>::SuperRef<double& (tempData::*)()> ref1(&tempData::getData);

    MemberReference<tempData, Vector<3,double>>::SuperRef ref2(std::make_pair(&tempData::getMomentum, &tempData::setMomentum));

    app = &ref;
    app = &ref1;
    tempData d;
    app->setValue(&d, 42.15);
    DBGVAR(app->getValue(&d));

    MemberApproach<tempData, Vector<3,double>>* app2 = &ref2;
    app2->setValue(&d, {42.15,84.30,42.15});
    DBGVAR(app2->getValue(&d), d.velocity);

    Traits<tempData, double, Vector<3,double>> r;
    Traits<tempData, double, Vector<3,double>>::getReference<0>() = &ref; //new MemberReference<tempData, double>::SuperRef<double tempData::*>(&tempData::density);
    Traits<tempData, double, Vector<3,double>>::getReference<1>() = new MemberReference<tempData, Vector<3, double>>::SuperRef(std::make_pair(&tempData::getMomentum, &tempData::setMomentum));
    DBGVAR(r.refs.MemRefs<0>::ref->getValue(&d), (Traits<tempData, double, Vector<3,double>>::getReference<1>()->getValue(&d)));
}

void testOrig() {
    Vertex<5, double> vert;
    vector<double> vec;
    DBGVAR(has_public_member<double>::value);
    DBGVAR(has_public_member<decltype (vert)>::value);
    DBGVAR(has_public_member<decltype(vec)>::value);
    DBGVAR(has_public_member<Temp>::value);

    Temp t;
    member_function_ptr<Temp, double, double>::type pom1 = &Temp::fun;
    member_ptr<Temp, double>::type pom = &Temp::data;
    DBGVAR(((&t)->*pom1)(0.0),pom);
    //auto op = &vector<double>::operator[];
    //member_const_function_ptr<vector<double>,const double&, size_t>::type c_at = &vector<double>::at;

    //member_function_ptr<vector<double>,double&, size_t>::type at = &vector<double>::at;

    //member_const_function_ptr<vector<double>,const double&, size_t>::type op = &vector<double>::operator[];

    //DBGVAR(is_same<decltype(&vector<double>::operator[]), typename member_const_function_ptr<vector<double>,const double&, size_t>::type>::value);

    DBGVAR(is_indexable<const vector<double>>::value);


    OutOfLineSpecialization<size_t>::printSeq(1);
    OutOfLineSpecialization<int>::printSeq(1);

    callSpec<int>();

}


template <typename T>
class Base{
public:
    T data;
    Base(T dat){
        data = dat;
    }
};

template <>
class Base<double>{
public:
    double data;
    Base(double dat){
        DBGMSG("double");
        data = dat;
    }
};



template <typename T1, typename T2>
class Base<std::pair<T1,T2>>{
public:
    T1 first;
    T2 second;
    Base(std::pair<T1,T2> pair){
        first = pair.first;
        second = pair.second;
    }
};


int main()
{

    Base b1(0.0);

    Base b2(std::pair<char,int>{'1',3});

    DBGVAR(b2.first,b2.second);
    testMemberRef();

    std::function<int(int)> fce = [&b1](int i){return b1.data + 42 + i;};

    std::cout << fce(2);

    return 0;
}
