// Test of Traits class
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
#include <list>
#include <map>
#include <array>
#include <string>
#include <type_traits>
#include "GTMesh/Traits/Traits.h"
#include <GTMesh/NumericStaticArray/Vector.h>

using double_tuple = std::tuple<double>;

MAKE_CUSTOM_TRAIT(
        std::tuple<double>,
        "double attr",
        std::make_pair(static_cast<const double&(*)(const std::tuple<double>&)>(std::get<0>),
                       static_cast<double&(*)(std::tuple<double>&)>(std::get<0>))
        );

struct testTuple{
    double attr1;
    int attr2;
    std::string attr3;
};

double& getData1(testTuple& s){return s.attr1;}

int& getData2(testTuple& s){return s.attr2;}

std::string& getData3(testTuple& s){return s.attr3;}


const double& getData1(const testTuple& s){return s.attr1;}

const int& getData2(const testTuple& s){return s.attr2;}

const std::string& getData3(const testTuple& s){return s.attr3;}


MAKE_CUSTOM_TRAIT_IO(
        testTuple,
        "1", std::make_pair(static_cast<const double&(*)(const testTuple&)>(getData1), static_cast<double&(*)(testTuple&)>(getData1)),
        "2", std::make_pair(static_cast<const int&(*)(const testTuple&)>(getData2), static_cast<int&(*)(testTuple&)>(getData2)),
        "3", std::make_pair(static_cast<const std::string&(*)(const testTuple&)>(getData3), static_cast<std::string&(*)(testTuple&)>(getData3))
        );

TEST( TupleTraitsTest, basicTest )
{
    EXPECT_FALSE((HasDefaultTraits<testTuple>::value));
    EXPECT_FALSE((HasDefaultArithmeticTraits<testTuple>::value));
    EXPECT_TRUE((HasDefaultIOTraits<testTuple>::value));
    testTuple tt;
    DefaultIOTraits<testTuple>::getTraits().setValue<1>(tt, 5);
    DefaultIOTraits<testTuple>::getTraits().setValue<2>(tt, "Hello :)");

    EXPECT_EQ( DefaultIOTraits<testTuple>::getTraits().getValue<1>(tt), 5 );
    EXPECT_EQ( DefaultIOTraits<testTuple>::getTraits().getValue<2>(tt), "Hello :)" );


    std::tuple<double> t{1.5};
    DefaultTraits<std::tuple<double>>::getTraits().getAttr<0>(t) = 2.5;
    EXPECT_EQ(DefaultTraits<std::tuple<double>>::getTraits().getAttr<0>(t), 2.5);
}


template< typename Real, size_t len >
class TemplateClass{
    std::array<Real, len> arr;

    friend class DefaultTraits<TemplateClass<Real, len>>;
};

MAKE_ATTRIBUTE_TEMPLATE_TRAIT( (TemplateClass<Real, len>),
                               (typename Real, size_t len),
                               arr );

TEST( TemplateTraitsTest, basicTest )
{
    TemplateClass<double, 3> tc;
    DefaultTraits<TemplateClass<double, 3>>::getTraits().getAttr<0>(tc) = {1.5,2.5,3.5};

    EXPECT_TRUE((HasDefaultTraits<TemplateClass<double, 3>>::value));
    EXPECT_TRUE((HasDefaultTraits<TemplateClass<int, 5>>::value));
    EXPECT_TRUE((IsDirectAccess<DefaultTraits<TemplateClass<int, 5>>::traitsType::memRefType<0>>::value));

    EXPECT_EQ((DefaultTraits<TemplateClass<int, 5>>::size()), 1);
    EXPECT_EQ((DefaultTraits<TemplateClass<int, 5>>::getTraits().getName<0>()), "arr");
    EXPECT_EQ((DefaultTraits<TemplateClass<double, 3>>::getTraits().getValue<0>(tc)), (std::array<double, 3>{1.5,2.5,3.5}));
}

template<typename Real>
class quantity {
public:
    Real rho;
    Real getPressure(Real T) const {
        return rho * T;
    }

    void setPressure(Real p, Real T) {
        rho = p / T;
    }
};

template <typename Real>
Real& Temperature() {
    static Real result =  5;
    return result;
}

template <typename Real>
Real getFunc(const quantity<Real>& obj){
    return obj.getPressure(Temperature<Real>());
}

template <typename Real>
void setFunc(quantity<Real>& obj, const Real& p){
    obj.setPressure(p, Temperature<Real>());
}

MAKE_CUSTOM_TEMPLATE_TRAIT((quantity<Real>), (typename Real),
                           "pressure",
                           std::make_pair( getFunc<Real>, setFunc<Real>)
                           );


template <typename Real>
auto bindPressure(Real& Temperature){
    return std::make_pair(
                std::bind(&quantity<Real>::getPressure, std::placeholders::_1, Temperature),
                std::bind(&quantity<Real>::setPressure, std::placeholders::_1, std::placeholders::_2, Temperature)
                );
}



TEST( TraitsOnBind, basicTest )
{

    quantity<double> instance;
    instance.rho = 1;
    EXPECT_EQ(DefaultTraits<quantity<double>>::getTraits().getValue<0>(instance), Temperature<double>());
    double T = 4 * Temperature<double>();
    // ability to bind std::bind using traits
    auto traits = makeTraits<quantity<double>>("4pressure", bindPressure<double>(T));

    EXPECT_EQ(traits.getValue<0>(instance), 4*Temperature<double>());
    EXPECT_TRUE((std::is_same<decltype(traits)::type<0>, double>::value));
    EXPECT_FALSE((IsDirectAccess<decltype(traits)::refType<0>>::value));

    // changing the value
    traits.setValue<0>(instance, 40);
    EXPECT_EQ(instance.rho, 2);
    EXPECT_EQ(DefaultTraits<quantity<double>>::getTraits().getValue<0>(instance), 10);
}

// Test binding lambda function
TEST( TraitsLambda, basicTest )
{

    quantity<double> instance;
    instance.rho = 1;

    double temperature = 50;
    auto traits = makeTraits<quantity<double>>(
        "pressure",
        std::make_pair(
            [&temperature](const quantity<double>& q){return q.getPressure(temperature);},
            [&temperature](quantity<double>& q, double p){return q.setPressure(temperature, p);}
        )
    );

    EXPECT_EQ(traits.getValue<0>(instance), 50);
    EXPECT_EQ(traits.getValue<0>(instance), instance.getPressure(temperature));
    EXPECT_TRUE((std::is_same<decltype(traits)::type<0>, double>::value));
    EXPECT_FALSE((IsDirectAccess<decltype(traits)::refType<0>>::value));


    // changing the value
    traits.setValue<0>(instance, 25);
    EXPECT_EQ(instance.rho, 2);
    EXPECT_EQ(DefaultTraits<quantity<double>>::getTraits().getValue<0>(instance), 10);
}
#endif

#include "UnitTests/main.h"
