#ifndef MACROFOREACH_H
#define MACROFOREACH_H

#define FOR_EACH_NARG(...) FOR_EACH_NARG_(__VA_ARGS__, FOR_EACH_RSEQ_N())
#define FOR_EACH_NARG_(...) FOR_EACH_ARG_N(__VA_ARGS__)
#define FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, N, ...) N
#define FOR_EACH_RSEQ_N() 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 09, 08, 07, 06, 05, 04, 03, 02, 01, 00


#define CONCATENATE(arg1, arg2)   CONCATENATE1(arg1, arg2)
#define CONCATENATE1(arg1, arg2)  CONCATENATE2(arg1, arg2)
#define CONCATENATE2(arg1, arg2)  arg1##arg2

#define FOR_EACH_00(what, ...)
#define FOR_EACH_01(what, x, ...)   what(x)
#define FOR_EACH_02(what, x, ...)   what(x), FOR_EACH_01(what, __VA_ARGS__)
#define FOR_EACH_03(what, x, ...)   what(x), FOR_EACH_02(what, __VA_ARGS__)
#define FOR_EACH_04(what, x, ...)   what(x), FOR_EACH_03(what, __VA_ARGS__)
#define FOR_EACH_05(what, x, ...)   what(x), FOR_EACH_04(what, __VA_ARGS__)
#define FOR_EACH_06(what, x, ...)   what(x), FOR_EACH_05(what, __VA_ARGS__)
#define FOR_EACH_07(what, x, ...)   what(x), FOR_EACH_06(what, __VA_ARGS__)
#define FOR_EACH_08(what, x, ...)   what(x), FOR_EACH_07(what, __VA_ARGS__)
#define FOR_EACH_09(what, x, ...)   what(x), FOR_EACH_08(what, __VA_ARGS__)
#define FOR_EACH_10(what, x, ...)   what(x), FOR_EACH_09(what, __VA_ARGS__)
#define FOR_EACH_11(what, x, ...)   what(x), FOR_EACH_10(what, __VA_ARGS__)
#define FOR_EACH_12(what, x, ...)   what(x), FOR_EACH_11(what, __VA_ARGS__)
#define FOR_EACH_13(what, x, ...)   what(x), FOR_EACH_12(what, __VA_ARGS__)
#define FOR_EACH_14(what, x, ...)   what(x), FOR_EACH_13(what, __VA_ARGS__)
#define FOR_EACH_15(what, x, ...)   what(x), FOR_EACH_14(what, __VA_ARGS__)
#define FOR_EACH_16(what, x, ...)   what(x), FOR_EACH_15(what, __VA_ARGS__)
#define FOR_EACH_17(what, x, ...)   what(x), FOR_EACH_16(what, __VA_ARGS__)
#define FOR_EACH_18(what, x, ...)   what(x), FOR_EACH_17(what, __VA_ARGS__)
#define FOR_EACH_19(what, x, ...)   what(x), FOR_EACH_18(what, __VA_ARGS__)
#define FOR_EACH_20(what, x, ...)   what(x), FOR_EACH_19(what, __VA_ARGS__)
#define FOR_EACH_21(what, x, ...)   what(x), FOR_EACH_20(what, __VA_ARGS__)
#define FOR_EACH_22(what, x, ...)   what(x), FOR_EACH_21(what, __VA_ARGS__)
#define FOR_EACH_23(what, x, ...)   what(x), FOR_EACH_22(what, __VA_ARGS__)
#define FOR_EACH_24(what, x, ...)   what(x), FOR_EACH_23(what, __VA_ARGS__)
#define FOR_EACH_25(what, x, ...)   what(x), FOR_EACH_24(what, __VA_ARGS__)
#define FOR_EACH_26(what, x, ...)   what(x), FOR_EACH_25(what, __VA_ARGS__)
#define FOR_EACH_27(what, x, ...)   what(x), FOR_EACH_26(what, __VA_ARGS__)
#define FOR_EACH_28(what, x, ...)   what(x), FOR_EACH_27(what, __VA_ARGS__)
#define FOR_EACH_29(what, x, ...)   what(x), FOR_EACH_28(what, __VA_ARGS__)
#define FOR_EACH_30(what, x, ...)   what(x), FOR_EACH_29(what, __VA_ARGS__)
#define FOR_EACH_31(what, x, ...)   what(x), FOR_EACH_30(what, __VA_ARGS__)
#define FOR_EACH_32(what, x, ...)   what(x), FOR_EACH_31(what, __VA_ARGS__)
#define FOR_EACH_33(what, x, ...)   what(x), FOR_EACH_32(what, __VA_ARGS__)
#define FOR_EACH_34(what, x, ...)   what(x), FOR_EACH_33(what, __VA_ARGS__)
#define FOR_EACH_35(what, x, ...)   what(x), FOR_EACH_34(what, __VA_ARGS__)
#define FOR_EACH_36(what, x, ...)   what(x), FOR_EACH_35(what, __VA_ARGS__)
#define FOR_EACH_37(what, x, ...)   what(x), FOR_EACH_36(what, __VA_ARGS__)
#define FOR_EACH_38(what, x, ...)   what(x), FOR_EACH_37(what, __VA_ARGS__)
#define FOR_EACH_39(what, x, ...)   what(x), FOR_EACH_38(what, __VA_ARGS__)
#define FOR_EACH_40(what, x, ...)   what(x), FOR_EACH_39(what, __VA_ARGS__)


#define FOR_EACH_(N, what, ...) CONCATENATE(FOR_EACH_, N)(what, __VA_ARGS__)
#define FOR_EACH(what, ...) FOR_EACH_(FOR_EACH_NARG(__VA_ARGS__), what, __VA_ARGS__)



#define FOR_EACH_EVEN_00(what, ...)
#define FOR_EACH_EVEN_02(what, x_odd, x_even, ...)   what(x_even)
#define FOR_EACH_EVEN_04(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_02(what, __VA_ARGS__)
#define FOR_EACH_EVEN_06(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_04(what, __VA_ARGS__)
#define FOR_EACH_EVEN_08(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_06(what, __VA_ARGS__)
#define FOR_EACH_EVEN_10(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_08(what, __VA_ARGS__)
#define FOR_EACH_EVEN_12(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_10(what, __VA_ARGS__)
#define FOR_EACH_EVEN_14(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_12(what, __VA_ARGS__)
#define FOR_EACH_EVEN_16(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_14(what, __VA_ARGS__)
#define FOR_EACH_EVEN_18(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_16(what, __VA_ARGS__)
#define FOR_EACH_EVEN_20(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_18(what, __VA_ARGS__)
#define FOR_EACH_EVEN_22(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_20(what, __VA_ARGS__)
#define FOR_EACH_EVEN_24(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_22(what, __VA_ARGS__)
#define FOR_EACH_EVEN_26(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_24(what, __VA_ARGS__)
#define FOR_EACH_EVEN_28(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_26(what, __VA_ARGS__)
#define FOR_EACH_EVEN_30(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_28(what, __VA_ARGS__)
#define FOR_EACH_EVEN_32(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_30(what, __VA_ARGS__)
#define FOR_EACH_EVEN_34(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_32(what, __VA_ARGS__)
#define FOR_EACH_EVEN_36(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_34(what, __VA_ARGS__)
#define FOR_EACH_EVEN_38(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_36(what, __VA_ARGS__)
#define FOR_EACH_EVEN_40(what, x_odd, x_even, ...)   what(x_even), FOR_EACH_EVEN_38(what, __VA_ARGS__)

#define FOR_EACH_EVEN_(N, what, ...) CONCATENATE(FOR_EACH_EVEN_, N)(what, __VA_ARGS__)
#define FOR_EACH_EVEN(what, ...) FOR_EACH_EVEN_(FOR_EACH_NARG(__VA_ARGS__), what, __VA_ARGS__)
/*
#define DUMMY
#define FOR_EACH_ODD_(N, what, ...) CONCATENATE(FOR_EACH_EVEN_, N)(what, __VA_ARGS__)
#define FOR_EACH_ODD(what, ...) FOR_EACH_EVEN_(FOR_EACH_NARG( , __VA_ARGS__), what, , __VA_ARGS__)
*/

#endif // MACROFOREACH_H
