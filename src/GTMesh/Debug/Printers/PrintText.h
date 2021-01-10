#ifndef PRINTTEXT_H
#define PRINTTEXT_H
#include <GTMesh/Traits/CustomTypeTraits.h>

struct PrintText {

    template<typename T>
    using isTextType = std::enable_if_t < std::is_same<T, std::string>::value ||
                                          std::is_same<T, const char*>::value ||
                                          std::is_same<T, char*>::value ||
                                          std::is_same<T, char>::value,
                                          bool >;

    template< typename TextType,
              isTextType<TextType> dummy = true >
    static void print(std::ostream& ost, const TextType& str, ...) {
        ost << '"' << str << '"';
    }

    static void print(const std::string& str, ...)
    {
        printf("\"%s\"", str.c_str());
    }


    static void print(const char* str, ...)
    {
        printf("\"%s\"", str);
    }


    static void print(const char str, ...)
    {
        printf("\"%c\"", str);
    }

};

#endif // PRINTTEXT_H
