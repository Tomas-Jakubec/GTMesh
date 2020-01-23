#define CONSOLE_COLORED_OUTPUT
#include "../src/Debug/Debug.h"
#include "../src/UnstructuredMesh/UnstructuredMesh.h"
#include "../src/Traits/MemberApproach/MemberApproach.h"
#include "../src/Traits/Traits.h"

#include <chrono>
#include <functional>
#include <type_traits>
#include <iostream>
#include <list>
#include <unordered_map>
#include <functional>
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
    double fun(double d [[maybe_unused]]){return data;}
};
#define STR(var) #var
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
    ConsoleLogger<>::writeVar(__LINE__, __FILE__, "r", r, "i", i, "c", c, "list", list, "vec", vec, "b", b, "map", m);
    ConsoleLogger<>::writeVar(__LINE__, __FILE__,"---", {5,4,3,2});
    DBGVAR(r, i, c, list, vec, b, m);

    Vertex<7, double> vert;
    DBGVAR(vert, vert.getCoordinates());

    DBGVAR((IsExportable<decltype(vec)>::value));

    DBGVAR((IsExportable<double>::value));

    DBGVAR(IsIndexable<double>::value);

    DBGVAR(IsIndexable<decltype(vec)>::value);

    DBGVAR(IsIndexable<decltype(vert)>::value);

    Subelement<size_t> s({1});
    auto v = {1,2,3};
    DBGVAR(s, v);

    DBGVAR_HTML(r, i, c, list, vec, b, m);

    DBGVAR_HTML(r+1, i+1, char(c+1), list, vec[0], b, m["prvni"]);
}



//The concept implementation
template<typename T>
class NeedIterator{
    static_assert (IsIterable<T>::value, "The type must be iterable");
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


template<typename T, typename Void = void>
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

    Vector<3,double> getMomentum()const{
        return velocity*density;
    }

    void setMomentum(const Vector<3,double>& val){
        velocity = val / density;
    }

};

/*
template<>
class Traits<tempData>{
public:
    using ttype = Traits<tempData, double, Vector<3,double>>;
    const static ttype tr;
};
const Traits<tempData>::ttype Traits<tempData>::tr("density", &tempData::density, "momentum"s, std::make_pair(&tempData::getMomentum, &tempData::setMomentum));
*/


//MAKE_NAMED_ATTRIBUTE_TRAIT(tempData, "density", density, "velocity", velocity);
//MAKE_ATTRIBUTE_TRAIT(tempData, density, velocity);

MAKE_CUSTOM_ATTRIBUTE_TRAIT(tempData, "density", &tempData::density, "momentum", std::make_pair(&tempData::getMomentum, &tempData::setMomentum))
MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(tempData, density, velocity);

struct ExportTest {
    int attrInt = 1;
    double attrDouble = 42.15;
    char attrChar = 42;
    std::string attrStr = "Ahojky";
    std::vector<std::string> attrVec = {"tohle", "je", "nejlepsi", "debugovaci", "system"};
    tempData attrTempData{42.15, {1,2,1}};
};
MAKE_ATTRIBUTE_TRAIT(ExportTest, attrInt, attrDouble, attrChar, attrStr, attrTempData, attrVec);
MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(ExportTest, attrInt, attrDouble, attrTempData);

#include "../src/Traits/TraitsAlgorithm/TraitsAlgorithm.h"

struct dataStruct {
    int iD;
    double dD;

}; MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(dataStruct, iD, dD);
void testMemberRef(){


    tempData d;

    //DBGVAR(Traits<tempData>::ttype::getName<0>());

    Traits<tempData>::getTraits().getReference<0>().setValue(&d, 0.0);
    DBGVAR(Traits<tempData>::getTraits().getReference<0>().getValue(&d));
    Traits<tempData>::getTraits().getReference<0>().setValue(d, 42.15);
    Traits<tempData>::getTraits().getReference<1>().setValue(&d, {42.15,84.30,42.15});

    DBGVAR(Traits<tempData>::getTraits().getName<0>(),(Traits<tempData>::getTraits().getReference<0>().getValue(&d)), Traits<tempData>::getTraits().getName<1>(), d.velocity);
    DBGVAR(Traits<tempData>::is_specialized,HasDefaultTraits<tempData>::value, d);

    ExportTest e;
    DBGVAR(e, ClassC<>());
}

void testTraitsAlgorithms() {
    ExportTest e1, e2;
    ExportTest res = e1 + e2;
    std::vector<ExportTest> vec(40, ExportTest());

    DBGVAR(2.45*e1,e1 + e2, e1, e2,HasDefaultArithmeticTraits<int>::value, max(e2), min(e2), max(vec));
}


struct NumStruct {
    double data;
};
MAKE_ATTRIBUTE_TRAIT(NumStruct, data);


struct NumStruct2 {
    double data1;
    double data2;
};
MAKE_ATTRIBUTE_TRAIT(NumStruct2, data1, data2);


template <unsigned int Index>
const MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<Index>, DefaultArithmeticTraits<tempData>::traitsType::refType<Index>>
getReference() {
    if constexpr (Index == 0){
        MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<0>, DefaultArithmeticTraits<tempData>::traitsType::refType<0>>
                mr1(&tempData::density);
        return mr1;
    }
    else{
        MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<1>, DefaultArithmeticTraits<tempData>::traitsType::refType<1>>
                mr2(&tempData::velocity);
        return mr2;
    }
}

struct References {
    template<unsigned int Index>
    using refType = typename DefaultArithmeticTraits<tempData>::traitsType::refType<Index>;

    template <unsigned int Index>
    using type = typename DefaultArithmeticTraits<tempData>::traitsType::type<Index>;


    template<unsigned int Index = 0, typename Dummy = void>
    struct MemRefs: public MemRefs<Index + 1> {

        const MemberReference<tempData, type<Index>, refType<Index>> ref;
        const char* name;


        template <typename ... REST>
        MemRefs(const char* n, refType<Index> r, REST... rest) : MemRefs<Index + 1> (rest...), ref(r), name(n){}
    };

    template<typename Dummy>
    struct MemRefs<1, Dummy>{
        const MemberReference<tempData, type<1>, refType<1>> ref;
        const char* name;

        MemRefs(const char* n, refType<1> r) : ref(r), name(n){}
    };

    const MemRefs<0, void> refs;

    References() : refs("",&tempData::density,"", &tempData::velocity){}

    template <unsigned int Index>
    const MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<Index>, DefaultArithmeticTraits<tempData>::traitsType::refType<Index>>
    getReference() {
        return refs.MemRefs<Index,void>::ref;
    }
};


struct traitPublisher {
    static const
    Traits<tempData, typename DefaultArithmeticTraits<tempData>::traitsType::refType<0>,
         typename DefaultArithmeticTraits<tempData>::traitsType::refType<1>> trait;

    static const
    Traits<tempData, typename DefaultArithmeticTraits<tempData>::traitsType::refType<0>,
    typename DefaultArithmeticTraits<tempData>::traitsType::refType<1>> getTrait(){
        return Traits<tempData, typename DefaultArithmeticTraits<tempData>::traitsType::refType<0>,
                typename DefaultArithmeticTraits<tempData>::traitsType::refType<1>>("1",&tempData::density,"2",&tempData::velocity);
    }

    static const
        MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<0>, DefaultArithmeticTraits<tempData>::traitsType::refType<0>>
                mr1;
    static const
        MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<1>, DefaultArithmeticTraits<tempData>::traitsType::refType<1>>
                mr2;
};
const
Traits<tempData, typename DefaultArithmeticTraits<tempData>::traitsType::refType<0>,
typename DefaultArithmeticTraits<tempData>::traitsType::refType<1>> traitPublisher::trait("1",&tempData::density,"2",&tempData::velocity);
const
        MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<0>, DefaultArithmeticTraits<tempData>::traitsType::refType<0>>
                traitPublisher::mr1(&tempData::density);
const
    MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<1>, DefaultArithmeticTraits<tempData>::traitsType::refType<1>>
            traitPublisher::mr2(&tempData::velocity);


void testNumericTraitsPerformance() {

    size_t size;
    DBGMSG("input size");
    cin >> size;

    std::vector<tempData> vec(size, {1, {1,1,1}});
    double ini;
    DBGMSG("input ini");
    std::cin >> ini;

    for (auto& val : vec) {
        val *= ini;
    }

    DBGMSG("input rep");
    int maxRep;
    cin >> maxRep;

    DBGVAR(size * ini * maxRep);



    auto clock = std::chrono::high_resolution_clock();
    DBGMSG("primary approach +=");
    long long deviation = 0;
    long long duration = 0;
    auto start = clock.now();
    tempData res({1,{1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res.density += vec[i].density;

            res.velocity += vec[i].velocity;

        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    auto avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("primary approach = +");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res.density = res.density + vec[i].density;

            res.velocity = res.velocity + vec[i].velocity;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += vec[i];
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DefaultArithmeticTraits<tempData>::traitsType::refType<0> const r1(&tempData::density);

    DefaultArithmeticTraits<tempData>::traitsType::refType<1> const r2(&tempData::velocity);

    DBGMSG("simulation numeric traits += (one ref one direct)");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            res.density += vec[i].*r1;

            res.velocity += vec[i].*r2;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    DBGMSG("simulation numeric traits += (both ref)");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            res.*r1 += vec[i].*r1;

            res.*r2 += vec[i].*r2;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("simulation numeric traits = +");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            res.*r1 = res.*r1 + vec[i].*r1;

            res.*r2 = res.*r2 + vec[i].*r2;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    DefaultArithmeticTraits<tempData>::traitsType::refType<0> volatile r1_1(&tempData::density);

    DefaultArithmeticTraits<tempData>::traitsType::refType<1> volatile r2_1(&tempData::velocity);


    DBGMSG("both volatile references = +");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            res.*r1 = res.*r1_1 + vec[i].*r1_1;

            res.*r2 = res.*r2 + vec[i].*r2_1;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("second reference is volatile = +");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            res.*r1 = res.*r1 + vec[i].*r1_1;

            res.*r2 = res.*r2 + vec[i].*r2_1;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<0>, DefaultArithmeticTraits<tempData>::traitsType::refType<0>>
            mr1(&tempData::density);
    MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<1>, DefaultArithmeticTraits<tempData>::traitsType::refType<1>>
            mr2(&tempData::velocity);


    DBGMSG("MemberReference +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            mr1.getAttr(res) += mr1.getAttr(vec[i]);

            mr2.getAttr(res) += mr2.getAttr(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("MemberReference = +");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            mr1.setValue(res, mr1.getValue(res) + mr1.getAttr(vec[i]));

            mr2.setValue(res, mr2.getValue(res) + mr2.getAttr(vec[i]));
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("published MemberReference +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            traitPublisher::mr1.getAttr(res) += traitPublisher::mr1.getAttr(vec[i]);

            traitPublisher::mr2.getAttr(res) += traitPublisher::mr2.getAttr(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    References rr;

    const MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<0>, DefaultArithmeticTraits<tempData>::traitsType::refType<0>>&
            rmr1 = rr.getReference<0>();//DefaultArithmeticTraits<tempData>::tr.getReference<0>();
    const MemberReference<tempData, DefaultArithmeticTraits<tempData>::traitsType::type<1>, DefaultArithmeticTraits<tempData>::traitsType::refType<1>>&
            rmr2 = rr.getReference<1>();//DefaultArithmeticTraits<tempData>::tr.getReference<1>();


    DBGMSG("reference to MemberReference +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            rmr1.getAttr(res) += rmr1.getAttr(vec[i]);

            rmr2.getAttr(res) += rmr2.getAttr(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("reference to MemberReference = +");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {

            rmr1.setValue(res, rmr1.getValue(res) + rmr1.getAttr(vec[i]));

            rmr2.setValue(res, rmr2.getValue(res) + rmr2.getAttr(vec[i]));
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;



    DBGMSG("numeric traits inplace direct approach of tr +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            DefaultArithmeticTraits<tempData>::getTraits().template getReference<0>().getAttr(res) +=
                    DefaultArithmeticTraits<tempData>::getTraits().template getReference<0>().getAttr(vec[i]);

            DefaultArithmeticTraits<tempData>::getTraits().template getReference<1>().getAttr(res) +=
                    DefaultArithmeticTraits<tempData>::getTraits().template getReference<1>().getAttr(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }
    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits inplace +=");

    const auto tr1 = DefaultArithmeticTraits<tempData>::getTraits();
    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            tr1.template getReference<0>().getAttr(res) +=
                    tr1.template getReference<0>().getAttr(vec[i]);

            tr1.template getReference<1>().getAttr(res) +=
                    tr1.template getReference<1>().getAttr(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    Traits<tempData, typename DefaultArithmeticTraits<tempData>::traitsType::refType<0>,
         typename DefaultArithmeticTraits<tempData>::traitsType::refType<1>>const trait("1",r1,"2",r2);

    DBGMSG("numeric traits inplace only RHS +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            res.density += trait.getReference<0>().getAttr(vec[i]);

            res.velocity += trait.getReference<1>().getAttr(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits inplace = +");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            res.density = trait.getReference<0>().getAttr(res) + trait.getReference<0>().getAttr(vec[i]);

            res.velocity = trait.getReference<1>().getAttr(res) + trait.getReference<1>().getAttr(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits inplace both referenced +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            trait.getReference<0>().getAttr(res) += trait.getReference<0>().getAttr(vec[i]);
            trait.getReference<1>().getAttr(res) += trait.getReference<1>().getAttr(vec[i]);

        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits static published both referenced +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            traitPublisher::trait.getReference<0>().getAttr(res) += traitPublisher::trait.getReference<0>().getAttr(vec[i]);
            traitPublisher::trait.getReference<1>().getAttr(res) += traitPublisher::trait.getReference<1>().getAttr(vec[i]);

        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits static published getTrait() both referenced +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            traitPublisher::getTrait().getReference<0>().getAttr(res) += traitPublisher::getTrait().getReference<0>().getAttr(vec[i]);
            traitPublisher::getTrait().getReference<1>().getAttr(res) += traitPublisher::getTrait().getReference<1>().getAttr(vec[i]);

        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    DBGMSG("numeric traits inplace +=");

    start = clock.now();
    res = tempData({1, {1,1,1}});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {


            trait.getReference<0>().setValue(res, res.density + trait.getReference<0>().getAttr(vec[i]));
            trait.getReference<1>().setValue(res, res.velocity + trait.getReference<1>().getAttr(vec[i]));


        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("simple numeric data test");

    vector<NumStruct> numVec(size,{ini});

    DBGMSG("direct approach +=");

    start = clock.now();
    NumStruct res1 = NumStruct({0});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < numVec.size(); i++) {

            res1.data += numVec[i].data;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res1, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits +=");

    start = clock.now();
    res1 = NumStruct({0});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < numVec.size(); i++) {

            res1 += numVec[i];
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res1, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;



    DBGMSG("simple numeric data test");

    vector<NumStruct2> numVec2(size,{ini,ini});

    DBGMSG("direct approach +=");

    start = clock.now();
    NumStruct2 res2 = NumStruct2({0,0});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < numVec2.size(); i++) {

            res2.data1 += numVec2[i].data1;

            res2.data2 += numVec2[i].data2;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res2, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;



    DBGMSG("Traits approach instead of numeric traits +=");
    Traits<NumStruct2> t;
    start = clock.now();
    res2 = NumStruct2({0,0});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < numVec2.size(); i++) {

            t.getTraits().getReference<0>().getAttr(res2) += t.getTraits().getReference<0>().getAttr(numVec2[i]);

            t.getTraits().getReference<1>().getAttr(res2) += t.getTraits().getReference<1>().getAttr(numVec2[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res2, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    DBGMSG("numeric traits +=");

    start = clock.now();
    res2 = NumStruct2({0,0});
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < numVec2.size(); i++) {

            res2 += numVec2[i];
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res2, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    DBGMSG("\n","ExportTest performance, direct approach");

    std::vector<ExportTest> vec3(size, ExportTest());

    start = clock.now();
    ExportTest res3 = ExportTest();
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec3.size(); i++) {

            res3.attrInt += vec3[i].attrInt;

            res3.attrDouble += vec3[i].attrDouble;

            res3.attrTempData.density = res3.attrTempData.density + vec3[i].attrTempData.density;

            res3.attrTempData.velocity = res3.attrTempData.velocity + vec3[i].attrTempData.velocity;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res3, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("numeric traits +=");

    start = clock.now();
    res3 = ExportTest();
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec3.size(); i++) {
            res3 += vec3[i];
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res3, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

}


/*
Test of trasposing vector of struct to struct of vectors
*/


template <typename DataType, typename DataTypeTrait = Traits<DataType>>
struct Container{
    template <unsigned int index = 0, typename Dummy = void>
    struct StructOfArrays : public StructOfArrays<index + 1>{
        std::vector<typename DataTypeTrait::template type<index>> vec;
    };


    template <typename Dummy>
    struct StructOfArrays<Traits<DataType>::traitsType::size() - 1, Dummy>{
        std::vector<typename DataTypeTrait::template type<DataTypeTrait::size() - 1>> vec;
    };

    StructOfArrays<> data;

    static constexpr unsigned int size() {
        return DataTypeTrait::size();
    }

    template <unsigned int pos>
    const char* name() {
        return Traits<DataType>::getTraits().template getName<pos>();
    }

    template <unsigned int pos>
    std::vector<typename Traits<DataType>::traitsType::template type<pos>>& getDataAtPos() {
        return data.StructOfArrays<pos>::vec;
    }
};


void testStructTransposition() {
    Container<ExportTest, Traits<ExportTest>::traitsType> data;

    data.getDataAtPos<5>().resize(2);
    DBGVAR(data.size(), data.getDataAtPos<5>(), data.name<5>());

};



/**
  Test TraitApply
  */
/*
template <typename T, typename Void = void>
struct TraitApply {

};
*/

/*
template < unsigned int index>
class TraitApply {
public:
    template <class Functor>
    static auto apply (Functor f,...)
    -> typename std::enable_if<std::is_assignable_v<
    std::function<
        void(unsigned int,
             std::unique_ptr<
                MemberApproach<
                    ExportTest,
                    Traits<ExportTest>::ttype::type<index>
                >>&,
             const std::string&
             )
    >, Functor>>::type
    {

        static_assert (std::is_assignable_v<
                std::function<
                    void(unsigned int,
                         std::unique_ptr<
                            MemberApproach<
                                ExportTest,
                                Traits<ExportTest>::ttype::type<index>
                            >>&,
                         const std::string&
                         )
                >, Functor>, "");

        f(index, Traits<ExportTest>::ttype::getReference<index>(), Traits<ExportTest>::ttype::getName<index>());
        TraitApply<index - 1>::apply(f);
    }



    template <class Functor>
    static auto apply (Functor f)
    -> typename std::enable_if<std::is_assignable_v<
    std::function<
        void(std::unique_ptr<
                MemberApproach<
                    ExportTest,
                    Traits<ExportTest>::ttype::type<index>
                >>&,
             const std::string&
             )
    >, Functor>>::type
    {

        static_assert (std::is_assignable_v<
                std::function<
                    void(std::unique_ptr<
                            MemberApproach<
                                ExportTest,
                                Traits<ExportTest>::ttype::type<index>
                            >>&,
                         const std::string&
                         )
                >, Functor>, "");

        f(Traits<ExportTest>::ttype::getReference<index>(), Traits<ExportTest>::ttype::getName<index>());
        TraitApply<index - 1>::apply(f);
    }

    template <template <typename, typename>class Functor>
    static auto apply ()
    -> typename std::enable_if<std::is_class<Functor<ExportTest, Traits<ExportTest>::ttype::type<index>>>::value>::type
    {

        static_assert (std::is_assignable<
                std::function<
                    void(unsigned int,
                         std::unique_ptr<
                            MemberApproach<
                                ExportTest,
                                Traits<ExportTest>::ttype::type<index>
                            >>&,
                         const std::string&
                         )
                >, Functor<ExportTest, Traits<ExportTest>::ttype::type<index>>>::value, "");


        Functor<ExportTest, Traits<ExportTest>::ttype::type<index>>()(index, Traits<ExportTest>::ttype::getReference<index>(), Traits<ExportTest>::ttype::getName<index>());
        TraitApply<index - 1>::template apply<Functor>();
    }

};


template <>
class TraitApply<0> {
public:
    template <class Functor>
    static auto apply (Functor f,...)
    -> typename std::enable_if<std::is_assignable_v<
    std::function<
        void(unsigned int,
             std::unique_ptr<
                MemberApproach<
                    ExportTest,
                    Traits<ExportTest>::ttype::type<0>
                >>&,
             const std::string&
             )
    >, Functor>>::type
    {
        f(0, Traits<ExportTest>::ttype::getReference<0>(), Traits<ExportTest>::ttype::getName<0>());
    }


    template <class Functor>
    static auto apply (Functor f)
    -> typename std::enable_if<std::is_assignable_v<
    std::function<
        void(std::unique_ptr<
                MemberApproach<
                    ExportTest,
                    Traits<ExportTest>::ttype::type<0>
                >>&,
             const std::string&
             )
    >, Functor>>::type
    {
        f(Traits<ExportTest>::ttype::getReference<0>(), Traits<ExportTest>::ttype::getName<0>());
    }

    template <template <typename, typename>class Functor>
    static auto apply ()
    -> typename std::enable_if<std::is_class<Functor<ExportTest,Traits<ExportTest>::ttype::type<0>>>::value>::type
    {
        Functor<ExportTest, Traits<ExportTest>::ttype::type<0>>()(0, Traits<ExportTest>::ttype::getReference<0>(), Traits<ExportTest>::ttype::getName<0>());
    }
};
*/

struct Depth {
    int value = 0;
};

template <typename Class, typename T>
class Func {
public:

    template <typename U = T>
    auto operator()(unsigned int index, const MemberApproach<Class, T>&, const std::string& name)
    -> typename std::enable_if<!(HasDefaultTraits<U>::value)>::type
    {
        DBGVAR(Singleton<Depth>::getInstance().value,index, name);
    }

    template <typename U = T>
    auto operator()(unsigned int index, const MemberApproach<Class, T>&, const std::string& name)
    -> typename std::enable_if<HasDefaultTraits<U>::value>::type
    {
        DBGVAR(Singleton<Depth>::getInstance().value,index, name);
        Singleton<Depth>::getInstance().value++;
        Traits<T>::getTraits().template apply<Func>();
        Singleton<Depth>::getInstance().value--;
    }
};



void testTraitApply() {

    auto lambda = [](unsigned int index, auto& , const std::string& name){DBGVAR(index, name);};

    DBGVAR(std::is_function<decltype (lambda)>::value);

    //TraitApply<5>::apply(lambda);

    auto lambda1 = []( auto& , const std::string& name){DBGVAR(name);};

    //TraitApply<5>::apply(lambda1);
DBGMSG("Tady");
    Traits<ExportTest>::getTraits().apply<Func>();
    //TraitApply<5>::apply<Func>();

    Traits<ExportTest>::getTraits().apply(lambda);

    Traits<ExportTest>::getTraits().apply(lambda1);

    //Traits<ExportTest>::ttype::apply<Func>();


}


/*
 *
 * Compile time traits DONE
 */

#include <cstdio>

struct Foo {
   int m;
   int r;

   int getM() {return m;}
   void setM(const int& _m) {m = _m;}
} foo = {2, 3};



template<typename R1, typename R2, R1 Mem1, R2 Mem2>
struct B {
   static constexpr R1 get = Mem1;

   static constexpr R2 set = Mem2;


   static typename MemberReferenceType<std::pair<R1,R2>>::type getValue(typename MemberReferenceType<std::pair<R1,R2>>::typeClass& foo) {
       return (foo.*(get))();
   }
};


template<typename Ref, Ref Mem>
struct B<Ref,Ref, Mem, Mem> {
   static constexpr Ref mp = Mem;

   inline static typename MemberReferenceType<Ref>::type& getValue(typename MemberReferenceType<Ref>::typeClass& foo) {
       return foo.*mp;
   }
};

#include<cstdio>
void testCompileTimeTraits() {
    typedef B<decltype(&Foo::m),decltype(&Foo::m), &Foo::m, &Foo::m> Bm;
    typedef B<decltype(&Foo::getM),decltype (&Foo::setM), &Foo::getM,&Foo::setM> Bp;
    typedef B<decltype(&Foo::r),decltype(&Foo::r), &Foo::r, &Foo::r> Br;
    DBGVAR(foo.*(Bm::mp), Bm::getValue(foo), foo.*(Br::mp), (foo.*(Bp::get))());
}

#include "../src/Singleton/Singleton.h"
/*
Test of order of constructors
*/
/*
struct mem{
    std::string s;
    mem(){DBGVAR(s);}
};

template<typename statMem>
class C1 {
    template<unsigned int Index>
    struct mem{
        std::string s;
        mem(){DBGVAR(s);}
    };

public:
    C1(){
        Singleton<mem<1>>::getInstance().s = "ahoj";
        DBGVAR("C1", Singleton<mem<1>>::getInstance().s);

    }
    static std::string& getS() {return Singleton<mem<1>>::getInstance().s;}

};


template <typename dummy>
class C2 {
public:
    static C1<mem> c;
    C2() {DBGVAR(c.getS());}
};
template <typename dummy> C1<mem> C2<dummy>::c;


template<typename statMem>
class C1_wrong {

public:

    static statMem s;

    C1_wrong(){
        s.s = "ahoj";
        DBGVAR("C1_wrong", s.s);

    }
    static std::string& getS() {return s.s;}

};
template <typename statMem> statMem C1_wrong<statMem>::s;



template <typename dummy>
class C2_wrong {
public:
    static C1_wrong<mem> c;
    C2_wrong() {DBGVAR(c.getS());}
};
template <typename dummy> C1_wrong<mem> C2_wrong<dummy>::c;



class C1_ {
public:
    static mem s1;
    C1_(){
        s1.s = "ahoj";
        DBGVAR("C1_", s1.s);

    }

};

mem C1_::s1;

class C2_ {
public:
    static C1_ c;
    C2_(){DBGVAR(c.s1.s);}
};
C1_ C2_::c;


using C = C1<mem>;

void testConstrucorOrder() {
    C2_ c;
    C2<void> c1;
    DBGVAR(C2_::c.s1.s, C2<void>::c.getS(), C::getS(), C2_wrong<void>::c.getS());
}
*/

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

class Number {
public:
    double num = 0;
    Number(double num) {
        this->num = num;
        DBGMSG("constructing number");
    }

    Number operator+(const Number& rhs){
        DBGMSG("operator+ const&");
        return Number(this->num + rhs.num);
    }


    Number operator+(Number&& rhs){
        DBGMSG("operator+ &&");
        rhs.num += this->num;
        return rhs;
    }

    Number operator*(const Number& rhs){
        return Number(this->num * rhs.num);
    }


    Number operator*(Number&& rhs){
        rhs.num *= this->num;
        return rhs;
    }


};

MAKE_ATTRIBUTE_TRAIT(Number,num);

void testOperator() {
    Number n(42.15), m(42);
    DBGMSG("start sum");
    Number nn = n + (n *(m + (m + n)));
    DBGVAR(nn);
}



template<typename Class, typename ValueType, typename Ref>
class TestMemberReference{

    TestMemberReference(Ref){}
};


template <typename Class, typename ValueType>
class TestMemberReference<Class, ValueType, ValueType Class::*> : MemberApproach<Class, ValueType>{

    using refType = ValueType Class::*;
public:
    refType const ref;

public:

    TestMemberReference(refType referenceToMember) : ref(referenceToMember){
        //ref = referenceToMember;
    }

    virtual ~TestMemberReference() = default;

    inline ValueType getValue(const Class* c) const {
        return c->*ref;
    }

    inline void setValue(Class* c, const ValueType& val) const {
        c->*ref = val;
    }

    inline ValueType getValue(const Class& c) const {
        return c.*ref;
    }

    inline void setValue(Class& c, const ValueType& val) const {
        c.*ref = val;
    }
};

template<typename Class, typename...RefTypes>
class TestTraits {
public:
    template<unsigned int Index>
    using refType = typename std::tuple_element<Index,std::tuple<RefTypes...>>::type;

    template <unsigned int Index>
    using type = typename MemberReferenceType<refType<Index>>::type;
private:
    template<unsigned int Index = 0, typename Dummy = void>
    struct MemRefs: public MemRefs<Index + 1> {

        const MemberReference<Class, type<Index>, refType<Index>> ref;
        std::string name;

        template <typename ... REST>
        MemRefs(std::string n, refType<Index> r, REST... rest) : MemRefs<Index + 1> (rest...), ref(r), name(n){}
    };

    template<typename Dummy>
    struct MemRefs<sizeof...(RefTypes) - 1, Dummy>{
        const MemberReference<Class, type<sizeof...(RefTypes) - 1>, refType<sizeof...(RefTypes) - 1>> ref;
        std::string name;

        MemRefs(std::string n, refType<sizeof...(RefTypes) - 1> r) : ref(r), name(n){}
    };

    const MemRefs<0, void> refs;
/*
    using refs = Singleton<MemRefs<sizeof... (RefTypes) - 1, void>>;

    template<unsigned int Pos = 0, typename ref, typename...Refs>
    static void _makeReferences(const std::string& name, ref member,Refs... refsAndNames) {
        _makeReferences<Pos, ref>(name, member);
        _makeReferences<Pos+1>(refsAndNames...);
    }

    template<unsigned int Pos, typename ref>
    static void _makeReferences(const std::string& name, ref member) {
        refs::getInstance().MemRefs<Pos, void>::name = name;
        refs::getInstance().MemRefs<Pos, void>::ref = MemberReference<Class, type<Pos>, refType<Pos>>(member);
    }


    template<unsigned int Pos, typename ref>
    static void _makeReferences(const char* name, ref member) {
        refs::getInstance().MemRefs<Pos, void>::name = name;
        refs::getInstance().MemRefs<Pos, void>::ref = MemberReference<Class, type<Pos>, refType<Pos>>(member);
    }
*/

public:



    static constexpr unsigned int size(){
        return sizeof... (RefTypes);
    }


    template<unsigned int Index>
    const MemberReference<Class, type<Index>, refType<Index>>& getReference(){
        return refs.MemRefs<Index, void>::ref;
    }

    template<unsigned int Index>
    type<Index> getValue(const Class* c){
        return getReference<Index>()->getValue(c);
    }

    template<unsigned int Index>
    type<Index> getValue(const Class& c){
        return getReference<Index>().getValue(c);
    }

    template<unsigned int Index>
    void setValue(Class* c, const type<Index>& val){
        getReference<Index>().setValue(c, val);
    }

    template<unsigned int Index>
    void setValue(Class& c, const type<Index>& val){
        getReference<Index>().setValue(c, val);
    }


    template<unsigned int Index>
    const std::string& getName(){
        return refs.MemRefs<Index, void>::name;
    }


    template<typename...Refs>
    TestTraits(Refs... refsAndNames) : refs(refsAndNames...){}


};


void testTestTraits(){

    TestTraits<ExportTest, decltype (&ExportTest::attrInt),  decltype (&ExportTest::attrDouble)> trait(
                "attrInt", &ExportTest::attrInt,
                "attrDouble", &ExportTest::attrDouble
                );
    ExportTest e;

    trait.setValue<0>(e,7);

    DBGVAR(trait.getValue<0>(e));

}


void testTraitPerformance() {

    size_t size;
    DBGMSG("input size");
    cin >> size;

    std::vector<ExportTest> vec(size);
    double ini;
    DBGMSG("input ini");
    std::cin >> ini;

    for (auto& val : vec) {
        val.attrDouble = ini;
    }

    DBGMSG("input rep");
    int maxRep;
    cin >> maxRep;

    DBGVAR(size * ini);

    auto clock = std::chrono::high_resolution_clock();
    DBGMSG("primary approach");
    long long deviation = 0;
    long long duration = 0;
    auto start = clock.now();
    double res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += vec[i].attrDouble;
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    auto avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;



    DBGMSG("constexpr ref");

    typedef B<decltype(&ExportTest::attrDouble),decltype(&ExportTest::attrDouble), &ExportTest::attrDouble, &ExportTest::attrDouble> doubleAttr;
    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += doubleAttr::getValue(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("direct ref");
    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += vec[i].*(doubleAttr::mp);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("member reference");
    //typedef B<decltype(&ExportTest::attrDouble),decltype(&ExportTest::attrDouble), &ExportTest::attrDouble, &ExportTest::attrDouble> doubleAttr;

    TestMemberReference<ExportTest, double, decltype (&ExportTest::attrDouble)> MR(&ExportTest::attrDouble);
    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += MR.getValue(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += vec[i].*(MR.ref);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("new test trait");
    TestTraits<ExportTest, decltype (&ExportTest::attrInt),  decltype (&ExportTest::attrDouble)> trait(
                "attrInt", &ExportTest::attrInt,
                "attrDouble", &ExportTest::attrDouble
                );

    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += trait.getValue<1>(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("member reference virtual");
    MemberReference<ExportTest, double, decltype (&ExportTest::attrDouble)> MR1(&ExportTest::attrDouble);
    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += MR1.getValue(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;

    DBGMSG("member reference virtual");
    const MemberApproach<ExportTest, double>* MA = new MemberReference<ExportTest, double, decltype (&ExportTest::attrDouble)>(&ExportTest::attrDouble);
    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += MA->getValue(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;


    DBGMSG("traits");
    start = clock.now();
    res = 0;
    for(int rep = 0; rep < maxRep; rep++){
        for(size_t i = 0; i < vec.size(); i++) {
            res += Traits<ExportTest>::getTraits().getValue<1>(vec[i]);
        }
        duration += (clock.now() - start).count();
        deviation += (clock.now() - start).count() * (clock.now() - start).count();
        start = clock.now();
    }

    avgDuration = duration / maxRep;
    DBGVAR(res, avgDuration , sqrt(((deviation)- maxRep * pow(avgDuration,2)) / (maxRep - 1)));
    deviation = 0;
    duration = 0;



}

void applyFunc(function_ptr<int, int>::type f, int arg) {
    DBGVAR(f(arg));
}
void applyFunc(std::function<int(int)>&f, int arg) {
    DBGVAR(f(arg));
}

template <typename Functor>
void applyFunc1(const Functor&f, int arg) {
    static_assert (std::is_assignable<std::function<void(int)>,Functor>::value,
                   "The Functor must be a function with one argument int and return type int");
    DBGVAR(f(arg));
}

void testFunction() {
    Base b1(0.0);

    Base b2(std::pair<char,int>{'1',3});

    DBGVAR(b2.first,b2.second);

    std::function<int(int)> fce ( [&b1](int i)->int{return b1.data + 42 + i;});

    applyFunc(fce, 2);
    applyFunc([](int i){return 42 + i;},2);
    applyFunc1([&b1](int i)->double{return b1.data + 42.15 + i;}, 2);
}



template <unsigned int dim, unsigned int Dim, ComputationMethod Method = DEFAULT>
struct calcCent{
    static void run() {DBGMSG("running default", dim);calcCent<dim + 1, Dim, Method>::run();}
};

template <unsigned int Dim, ComputationMethod Method>
struct calcCent<Dim, Dim, Method>{
    static void run() {DBGMSG("running default", Dim);}
};

template <unsigned int Dim, ComputationMethod Method>
struct calcCent<0, Dim, Method>{
    static void run() {DBGMSG("running default", 0);calcCent<1, Dim, Method>::run();}
};

template <unsigned int Dim, ComputationMethod Method>
struct calcCent<1, Dim, Method>{
    static void run() {DBGMSG("running default", 1);calcCent<2, Dim, Method>::run();}
};

template <>
struct calcCent<2, 3, TESSELLATED>{
    static void run() {DBGMSG("running MESH_TESSELLATED");calcCent<3, 3, TESSELLATED>::run();}
};


void testCalcCent() {
    calcCent<0,3>::run();
    calcCent<0,3, TESSELLATED>::run();
}


/*
    test of custom hash
*/
namespace std {
template<typename IndexType>
class hash<std::vector<IndexType>>{
public:
    size_t operator()(const std::vector<IndexType>& v) const {

        std::string_view sv(reinterpret_cast<const char*>(v.data()), sizeof(IndexType) * v.size());
        return std::hash<std::string_view>{}(sv);
    }
};

template<unsigned int Dim, typename Real>
class hash<Vertex<Dim, Real>>{
public:
    size_t operator()(const Vertex<Dim, Real>& vert) const {

        std::string_view sv(reinterpret_cast<const char*>(&vert), sizeof(Real) * vert.size());
        return std::hash<std::string_view>{}(sv);
    }

};
}

void testCustomUnorderedMap() {
    std::unordered_map<std::vector<size_t>, size_t> m;

    std::vector<size_t> v = {1,2,3};
    std::vector<size_t> v2 = {3,2,1};

    std::sort(v2.begin(), v2.end());
    DBGVAR(v == v2, sizeof(std::vector<size_t>), sizeof(std::string));


    DBGVAR(m[v] = 1, m[v2]);

    Vertex<3, double> vert1 = {1,2,3};

    std::unordered_map<Vertex<3,double>, size_t> mm;

    DBGVAR(mm[vert1] = 3);
    DBGVAR(*mm.begin(), sizeof(decltype(vert1)));


    std::set<int> s = {1,15,6,8};

    std::vector<int> vec(4);
    std::vector<int> vec2;

    std::copy(s.begin(), s.end(), vec.begin());

    //std::copy(s.begin(), s.end(), std::inserter(vec2, vec2.begin()));

    vec2.insert(vec2.begin(),s.begin(), s.end());

    std::vector<int> vec3(s.begin(), s.end());

    DBGVAR(vec,s, vec2, vec.capacity(), vec2.capacity(), vec3);
}


struct privateAttr{
private:
    std::string attr = "attr";
public:
    const std::string& getAttr(){return attr;}
    friend Traits<privateAttr>;
};

MAKE_NAMED_ATTRIBUTE_TRAIT(privateAttr, "str_attr", attr);


void testPrivateTrait(){
    privateAttr a;
    DBGVAR(a);
    Traits<privateAttr>::getTraits().setValue<0>(a, "new value");
    DBGVAR(a);

}

/*
#include "json.hpp"


struct person{
   std::string name, surname;

   struct address {
       std::string address;
       int num;
   } addr;
};

MAKE_ATTRIBUTE_TRAIT(person::address, address, num);
MAKE_ATTRIBUTE_TRAIT(person, name, surname, addr);

using json = nlohmann::json;

template <unsigned int Index>
struct TraitToJson {
    template<typename traitedClass>
    static
    typename enable_if<
        (Traits<traitedClass>::ttype::size() - 1 > Index)
    >::type
    to_json(json& j, const traitedClass& t) {

//        j.at(Traits<traitedClass>::ttype::template getName<Index>()) = Traits<traitedClass>::ttype::template getValue<Index>(t);
        j.push_back({Traits<traitedClass>::ttype::template getName<Index>(), Traits<traitedClass>::ttype::template getValue<Index>(t)});
        TraitToJson<Index + 1>::to_json(j, t);
    }


    template<typename traitedClass>
    static
    typename enable_if<
        (Traits<traitedClass>::ttype::size() - 1 == Index)
    >::type
    to_json(json& j, const traitedClass& t) {

        j.push_back({Traits<traitedClass>::ttype::template getName<Index>(), Traits<traitedClass>::ttype::template getValue<Index>(t)});
        //j.at(Traits<traitedClass>::ttype::template getName<Index>()) = Traits<traitedClass>::ttype::template getValue<Index>(t);

    }
};

template<typename traitedClass>
typename
std::enable_if<
    HasDefaultTraits<traitedClass>::value
>::type
to_json(json& j, const traitedClass& t){
    j = json::object();
    TraitToJson<0>::to_json(j, t);
}


template <unsigned int Index>
struct JsonToTrait {
    template<typename traitedClass>
    static
    typename enable_if<
        (Traits<traitedClass>::ttype::size() - 1 > Index)
    >::type
    from_json(const json& j, traitedClass& t) {

        Traits<traitedClass>::ttype::template setValue<Index>(
            t,
            j.at(Traits<traitedClass>::ttype::template getName<Index>()).template get<typename Traits<traitedClass>::ttype::template type<Index>>()
        );
        JsonToTrait<Index + 1>::from_json(j, t);
    }


    template<typename traitedClass>
    static
    typename enable_if<
        (Traits<traitedClass>::ttype::size() - 1 == Index)
    >::type
    from_json(const json& j, traitedClass& t) {

        Traits<traitedClass>::ttype::template setValue<Index>(
            t,
            j.at(Traits<traitedClass>::ttype::template getName<Index>()).template get<typename Traits<traitedClass>::ttype::template type<Index>>()
        );

    }
};

template<typename traitedClass>
typename
std::enable_if<
    HasDefaultTraits<traitedClass>::value
>::type
from_json(const json& j, traitedClass& t) {
    JsonToTrait<0>::from_json(j,t);
}

namespace ns {
    // a simple struct to model a person
    struct person {
        std::string name;
        std::string address;
        int age;
    };
}
namespace ns {
    void to_json(json& j, const person& p) {
        j = json{{"name", p.name}, {"address", p.address}, {"age", p.age}};
    }

    void from_json(const json& j, person& p) {
        j.at("name").get_to(p.name);
        j.at("address").get_to(p.address);
        j.at("age").get_to(p.age);
    }
} // namespace ns


void testJson() {

    auto js = R"({"name":"Tomik","surname":"Jakubec"})"_json;


    // create a person
    ns::person p {"Ned Flanders", "744 Evergreen Terrace", 60};

    // conversion: person -> json
    json j = p;

    std::cout << j << std::endl;
    // {"address":"744 Evergreen Terrace","age":60,"name":"Ned Flanders"}

    // conversion: json -> person
    auto p2 = j.get<ns::person>();

    person p_test = {"tomik","...", {"MS", 334}};

    json j_test = p_test;


    j_test.at("name") = "Ivisek";

    p_test = j_test;

    std::cout << j_test << std::endl;
    DBGVAR(j, p_test, j_test);

}*/






int main()
{
    //testDebug();
    //testOperator();
    //testMemberRef();
    //testConstrucorOrder();
    //testFunction();
    //testCalcCent();
    //testStructTransposition();
    //testTraitApply();
    //testCompileTimeTraits();
    //testTraitPerformance();
    //testCustomUnorderedMap();
    //testPrivateTrait();
    //testJson();
    //testTestTraits();
    //testTraitsAlgorithms();
    testNumericTraitsPerformance();
    return 0;
}


/** GCC error
#include <iostream>
#include <string>

struct mem{
    std::string s;
    mem(){std::cout <<"mem constructor: " << s << std::endl;}
};

template<typename statMem>
class C1_wrong {

public:

    static statMem s;

    C1_wrong(){
        s.s = "hello";
        std::cout << "C1_wrong s.s: " << s.s << std::endl;

    }
    static std::string& getS() {return s.s;}

};
template <typename statMem> statMem C1_wrong<statMem>::s;



template <typename dummy>
class C2_wrong {
public:
    static C1_wrong<mem> c;
    C2_wrong() {std::cout << "C2_wrong constructor c.getS(): " << c.getS() << std::endl;}
};
template <typename dummy> C1_wrong<mem> C2_wrong<dummy>::c;


void testConstrucorOrder() {

    std::cout << "global value of C2: " << C2_wrong<void>::c.getS() << std::endl;
}

int main()
{
    testConstrucorOrder();
}
  */
