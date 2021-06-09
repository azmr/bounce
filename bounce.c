// TODO: is there a way to do default args by this approach? (i.e. if the last arg is not provided it gets set to default)
// TODO: add a fixup arg type, where the instruction memory of the value be set is returned so that it can be altered by the caller
//          - this could default to 0 (don't need a larger input array) or could include passing a default value
//          - bnc_copy could actually use this by only passing a size and having the caller fill out the data
//          - careful with thread safety etc
//          - or could overwrite all array values with their eventual location & size
// TODO: type checking
// TODO: manual tail-call optimization? https://ekins.space/2017/06/07/trampoline/

typedef enum BncReg {
    Bnc_ax,     // rax, eax, ax
    Bnc_cx,
    Bnc_dx,
    Bnc_r8 = 8, // r8,  r8d, r8w
    Bnc_r9,
    Bnc_r10,
    Bnc_r11, // NOTE: volatile & tmp_reg for both Windows & Linux => good scratch register

    Bnc_xmm0 = 0,
    Bnc_xmm1,
    Bnc_xmm2,
    Bnc_xmm3,
    Bnc_xmm4,
    Bnc_xmm5,
    Bnc_xmm6,
    Bnc_xmm7,
} BncReg;

typedef enum BncTag {
    Bnc_fix       = (1 << 31),
    Bnc_float     = (1 << 30),
    Bnc_copy      = (1 << 29),
    Bnc_size_mask = Bnc_copy - 1,
    Bnc_reg       = 0,
} BncTag;

typedef union  { void *ptr; U8 rU8; U4 rU4; S8 rS8; F4 xF4; F8 xF8; BncReg reg; } BncVal;
typedef struct { BncTag tag; union{ BncVal; BncVal val; }; } BncArg; // tag contains size

internal inline BncArg bnc_copy(void *v, U4 size)
{
    assert(size <= Bnc_size_mask, "size must be less than %zu. size = %zu", Bnc_size_mask, size);
    return (BncArg){ Bnc_copy | Bnc_fix | size, .ptr = v };
}
internal inline BncArg bnc_struct(void *v, U4 size)
{
    BncArg result = {0};
    if (size > sizeof(void*)) { result = bnc_copy(v, size);                                  }
    else                      { result.tag = (Bnc_fix | size), memcpy(&result.rU8, v, size); }
    return result;
}

internal inline BncArg bnc_ptr(void *v) { return (BncArg){ (BncTag)(Bnc_fix | sizeof(void *)),         .ptr = v }; }
internal inline BncArg bnc_U8(U8     v) { return (BncArg){ (BncTag)(Bnc_fix | sizeof(U8)),             .rU8 = v }; }
internal inline BncArg bnc_S8(S8     v) { return (BncArg){ (BncTag)(Bnc_fix | sizeof(S8)),             .rS8 = v }; }
internal inline BncArg bnc_F4(F4     v) { return (BncArg){ (BncTag)(Bnc_float | Bnc_fix | sizeof(F4)), .xF4 = v }; }
internal inline BncArg bnc_F8(F8     v) { return (BncArg){ (BncTag)(Bnc_float | Bnc_fix | sizeof(F8)), .xF8 = v }; }


#if 1  // MACHINE CODE
internal inline Size mc_cpy(Byte **cur, Byte const *mc, Size size)
{
    if (*cur)
    {   memcpy(*cur, mc, size), *cur += size;   }
    return size;
}

#if 1  // MOV
internal inline Size mov_reg_lit4(Byte **cur, BncReg reg, BncVal src, BncReg tmp_reg)
{   (void)tmp_reg;
    U4   v       = src.rU4;
    Size mc_size = ((reg >= Bnc_r8) ? mc_cpy(cur, &(Byte){ 0x41 }, 1) : 0); // extended prefix

    Byte mc[] = {
        0xb8 + (reg % Bnc_r8),                        // mov literal into r
        (v&255), (v>>8)&255, (v>>16)&255, (v>>24)&255, // LE 4-byte value
    };
    mc_size += mc_cpy(cur, mc, sizeof(mc));

    return mc_size;
}

internal inline Size mov_reg_imm(Byte **cur, BncReg reg, BncVal src, BncReg tmp_reg)
{   (void)tmp_reg;
    U8   v    = src.rU8;
    Byte mc[] = {
        0x48 + (reg >= Bnc_r8),                            // 8-byte modifier prefix
        0xb8 + (reg %  Bnc_r8),                            // mov literal into r
        (v&255),     (v>>8)&255,  (v>>16)&255, (v>>24)&255, // LE 8-byte value
        (v>>32)&255, (v>>40)&255, (v>>48)&255, (v>>56)&255,
    };
    return mc_cpy(cur, mc, sizeof(mc));
}


internal inline Size mov_reg_reg(Byte **cur, BncReg dst, BncVal src, BncReg tmp_reg)
{   (void)tmp_reg;
    if (dst != src.reg)
    {
        // NOTE: 4-byte version takes 2 bytes instead of 3 for e_x but not r8d,r9d
        // given that we can only swap these low registers max once per function, it's not really worth the complexity for 1 byte...
        Byte mc[] = {
            (0x48 + (src.rU4 >= Bnc_r8) + 4*(dst >= Bnc_r8)), // 8-byte mod prefix
            0x89,                                                   // mov from register
            (0xc0 + 8*(src.rU4 % 8) + (dst % 8)),                   // src & dst reg
        };
        return mc_cpy(cur, mc, sizeof(mc));
    }
    else return 0;
}

internal inline Size mov_xmm_xmm(Byte **cur, BncReg dst_xmm, BncVal src_xmm, BncReg tmp_reg)
{   (void)tmp_reg;
    if (dst_xmm != src_xmm.reg)
    {
        Byte mc[] = {
            0xf3, 0x0f,                   // prefix
            0x7e,                         // mov to xmm
            (0xc0 + 8*dst_xmm + src_xmm.rU4), // src & dst
        };
        return mc_cpy(cur, mc, sizeof(mc));
    }
    else return 0;
}

// TODO: there may be a better approach; this is quite a large instruction
internal inline Size mov_xmm_F4(Byte **cur, int xmm, BncVal val, BncReg tmp_reg)
{
    Size mc_size = mov_reg_lit4(cur, Bnc_ax, val, tmp_reg);
    Byte mc[]    = {
        0x66, 0x0f,       // prefix
        0x6e,             // mov to xmm
        (0xc0 + 8*(xmm)), // dst reg
    };
    return mc_size + mc_cpy(cur, mc, sizeof(mc));
}

// TODO: there may be a better approach; this is quite a large instruction
internal inline Size mov_xmm_F8(Byte **cur, int xmm, BncVal val, BncReg tmp_reg)
{   (void)tmp_reg;
    assert(tmp_reg == Bnc_ax);
    Size mc_size = mov_reg_imm(cur, Bnc_ax, val, tmp_reg);
    Byte mc[]    = {
        0x66, 0x48, 0x0f, // prefix
        0x6e,             // mov to xmm
        (0xc0 + 8*(xmm)), // dst reg
    };
    return mc_size + mc_cpy(cur, mc, sizeof(mc));
}

typedef enum { Bnc_reg_mem, Bnc_mem_reg } BncMemRegDir;
internal inline Size mov_reg_mem_reg(Byte **cur, BncMemRegDir dir, BncReg reg, Size rsp_offset)
{
    Byte dir_mc[] = { [Bnc_reg_mem]=0x8b, [Bnc_mem_reg]=0x89, };
    if (0 <= rsp_offset && rsp_offset < 0x7f)
    {
        Byte mc[] = {
            0x48 + 4*(reg >= Bnc_r8),                   // prefix
            dir_mc[dir], 0x44 + 8*(reg % Bnc_r8), 0x24, // mov to reg
            (S1)rsp_offset                               // from [rsp + rsp_offset]
        }; // add rsp -size
        return mc_cpy(cur, mc, sizeof(mc));
    }
    else
    {
        TODO("handle larger rsp offsets (%zu)", rsp_offset);
        return 0;
    }
}
internal inline Size mov_reg_mem(Byte **cur, BncReg reg,            BncVal rsp_offset,     BncReg tmp_reg) { (void)tmp_reg; return mov_reg_mem_reg(cur, Bnc_reg_mem, reg,     rsp_offset.reg); }
internal inline Size mov_mem_reg(Byte **cur, int    rsp_offset,     BncVal reg,            BncReg tmp_reg) { (void)tmp_reg; return mov_reg_mem_reg(cur, Bnc_mem_reg, reg.reg, rsp_offset); }
internal inline Size mov_mem_mem(Byte **cur, int    dst_rsp_offset, BncVal src_rsp_offset, BncReg tmp_reg)
{   (void)tmp_reg;
    Size mc_size  = 0;
    if (dst_rsp_offset != src_rsp_offset.reg)
    {
        mc_size += mov_reg_mem_reg(cur, Bnc_reg_mem, tmp_reg, src_rsp_offset.reg);
        mc_size += mov_reg_mem_reg(cur, Bnc_mem_reg, tmp_reg, dst_rsp_offset);
    }
    return mc_size;
}

internal inline Size mov_mem_xmm(Byte **cur, int rsp_offset, BncVal xmm, BncReg tmp_reg)
{   (void)tmp_reg;
    if (0 <= rsp_offset && rsp_offset < 0x7f)
    {
        Byte mc[] = {
            0x66, 0x0f,             // prefix
            0xd6, 0x44 + 8*xmm.reg, // mov from xmm
            0x24,
            (S1)rsp_offset          // to [rsp + rsp_offset]
        }; // add rsp -size
        return mc_cpy(cur, mc, sizeof(mc));
    }
    else
    {
        TODO("handle larger rsp offsets (%zu)", rsp_offset);
        return 0;
    }
}


internal inline Size mov_mem_imm(Byte **cur, int rsp_offset, BncVal val, BncReg tmp_reg)
{
    // TODO (perf,mem): can do much more succinctly with values that fit in 4 bytes
    Size   mc_size = mov_reg_imm(cur, tmp_reg,    val,                    tmp_reg);
    return mc_size + mov_mem_reg(cur, rsp_offset, (BncVal){.reg=tmp_reg}, tmp_reg);
}

#endif // MOV

internal inline Size sub_rsp(Byte **cur, Size size)
{
    // NOTE: we negate this to get the largest negative 1-byte representation and add it;
    // we add a negative value rather than sub a positive one because
    // this way we get to pass 1 more parameter in 1 byte
    if (size < 0x80)
    {
        Byte mc[] = { 0x48, 0x83, 0xc4, (Byte)-(S1)size }; // add rsp -size
        return mc_cpy(cur, mc, sizeof(mc));
    }

    else TODO("haven't handled multi-byte sizes; trying to pass in %zu args", size / 8);
    return 0;
}

internal inline Size add_rsp(Byte **cur, Size size)
{
    // NOTE: we negate this to get the largest negative 1-byte representation and add it;
    // we sub a negative value rather than add a positive one because
    // this way we get to pass 1 more parameter in 1 byte
    if (size <= 0x80)
    {
        Byte mc[] = { 0x48, 0x83, 0xec, (Byte)-(S1)size }; // sub rsp -size
        return mc_cpy(cur, mc, sizeof(mc));
    }

    else TODO("haven't handled multi-byte sizes; trying to pass in %zu args", size / 8);
    return 0;
}

internal inline Size jmp_fn(Byte **cur, void *fn)
{
    Size mc_size = mov_reg_imm(cur, Bnc_ax, (BncVal){.ptr=fn}, 0); // mov rax, fn
    if (cur)
    {   cur[0][0]=0xff, cur[0][1]=0xe0; *cur += 2;   } // jmp rax
    return mc_size + 2;
}

internal inline Size call_fn(Byte **cur, void *fn)
{
    Size mc_size = mov_reg_imm(cur, Bnc_ax, (BncVal){.ptr=fn}, 0); // mov rax, fn
    if (cur)
    {   cur[0][0]=0xff, cur[0][1]=0xd0; *cur += 2;   } // call rax
    return mc_size + 2;
}

internal inline Size ret(Byte **cur)
{
    if (*cur) { **cur = 0xc3; ++*cur; } // ret
    return 1;
}

#endif // MACHINE CODE

static inline U1
spec_ceil_pow2_U1(U1 v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v++;
    return v;
}
internal inline UPtr spec_align_up_offset(Byte *base, UPtr offset, UPtr size)
{
    Size align = size >= 16 ? 16 : spec_ceil_pow2_U1((U1)size);
    assert(is_pow_2(align));
    --align;
    UPtr new_base = (UPtr)(base + offset + align) & ~align;
    return new_base - (UPtr)base;
}


#if _WIN32
// TODO: structs > 64 bits have a pointer passed to them in rcx
// if exe is null, cast the result to UPtr to get the size needed // TODO: or
// NOTE: the number & ordering of args corresponds to those in the target (not generated) function
// TODO: rename? trampoline? bounce? spring? preset_fn_args
internal void *
specialize(Byte *mem, void *fn, BncArg const args[], Size args_n)
{
#define assume_args_are_sorted
    persist const BncReg Arg_I_Regs[] = { Bnc_cx, Bnc_dx, Bnc_r8, Bnc_r9, };
    typedef Size BncMovFn(Byte **, int, BncVal, int);

    enum         { Loc_null, Loc_reg,     Loc_xmm,     Loc_mem,                                              Loc_Count };
    typedef enum { Dst_null, Dst_reg,     Dst_xmm,     Dst_mem,                                              Dst_Count } BncDstLoc;
    typedef enum { Src_null, Src_reg,     Src_xmm,     Src_mem,     Src__U8,     Src__F4,     Src__F8,       Src_Count } BncSrcLoc;
    persist BncMovFn * const Movs[Dst_Count][Src_Count] = {
        [Dst_reg]={ 0,       mov_reg_reg, 0,           0,           mov_reg_imm, },
        [Dst_xmm]={ 0,       0,           mov_xmm_xmm, 0,           0,           mov_xmm_F4,  mov_xmm_F8  },
        [Dst_mem]={ 0,       mov_mem_reg, mov_mem_xmm, mov_mem_mem, mov_mem_imm, mov_mem_imm, mov_mem_imm },
    };
    persist int const Arg_Locs[2/*stack*/][2/*float*/] = { // NOTE: other than immediates
        //         scalar   float
        /* regs */ Loc_reg, Loc_xmm,
        /* stack*/ Loc_mem, Loc_mem,
    };
    persist BncSrcLoc const Imm_Locs[2/*stack*/][2/*size*/][2/*float*/] = {
        //          scalar   float
        { // regs
            /* 4 */ Src__U8, Src__F4,
            /* 8 */ Src__U8, Src__F8, // only worth dealing with floats specially when they're going into XMM registers
        },
        { // stack
            /* 4 */ Src__U8, Src__U8,
            /* 8 */ Src__U8, Src__U8,
        }
    };

    Byte *data      = mem;
    Size  mc_size   = 0,
          data_size = 0,
          data_used = 0;

    Size gen_fn_args_n = 0;
    for (Size i = args_n; i-- > 0;)
    {
        BncArg arg     = args[i];
        gen_fn_args_n += !(arg.tag & Bnc_fix); // count the number of args to be passed to the generated fn
        if (arg.tag & Bnc_copy)
        { // determine the space needed before the code for copied data
            U4 arg_size = (arg.tag & Bnc_size_mask);
            data_size   = spec_align_up_offset(data, data_size, arg_size) + arg_size;
        }
    }
    Size src_fn_arg_i = gen_fn_args_n - 1;

    Byte *exe = (mem ? mem + data_size : 0); // make space for data before the executable
    Byte *cur = exe;

    Size rsp_d = 0;
    if (args_n > 4)
    {
        rsp_d    = args_n - gen_fn_args_n;
        rsp_d    = (rsp_d + rsp_d & 1) * sizeof(void*); // rsp 16 byte align
        mc_size += sub_rsp(&cur, rsp_d);
        mc_size += mov_mem_mem(&cur, 0, (BncVal){.reg=rsp_d}, Bnc_ax); // move the return address to rsp
    }


    // NOTE: need to make sure we don't clobber args before we use them elsewhere
    for (Size dst_fn_arg_i = args_n; dst_fn_arg_i-- > 0;)
    { // move args from a literal/generated position to their target position
        BncArg arg      = args[dst_fn_arg_i];
        U4     arg_size = (arg.tag & Bnc_size_mask) ?: 8; // NOTE: default 0 allows API users to not specify all args that haven't changed (as long as they're scalar)
        int    is_float = !!(arg.tag & Bnc_float);

        // COPY CACHED DATA ////////////////////////////////////////////////////////////
        if (arg.tag & Bnc_copy) // NOTE: this has to be done before using the arg's val
        { // push some data and pass a pointer to it
            data_used = spec_align_up_offset(data, data_used, arg_size);
            if (cur)
            {   arg.ptr = memcpy(&data[data_used], arg.ptr, arg_size);   }
            data_used += arg_size;
            arg_size  = sizeof(void*);
        }

        // DETERMINE PARAMETERS BASED ON ARG LOCATION //////////////////////////////////
        U4 const dsts[Dst_Count] = {
            [Dst_reg] = Arg_I_Regs[dst_fn_arg_i],
            [Dst_xmm] = dst_fn_arg_i,
            [Dst_mem] = 8*dst_fn_arg_i,
        };
        BncVal const srcs[Src_Count] = {
            [Src_reg] = { .rU4 = Arg_I_Regs[src_fn_arg_i] },
            [Src_xmm] = { .rU4 = src_fn_arg_i },
            [Src_mem] = { .rU4 = 8*src_fn_arg_i },
            [Src__U8] = arg.val,
            [Src__F4] = arg.val,
            [Src__F8] = arg.val,
        };

        // DETERMINE DEST & SRC LOCATIONS //////////////////////////////////////////////
        BncDstLoc dst = Arg_Locs[dst_fn_arg_i >= 4][is_float];
        BncSrcLoc src = 0;

        if (arg.tag & Bnc_fix)
        {
            assert(arg_size <= 8, "unhandled arg size: %u", arg_size);
            assert(!is_float || arg_size == 4 || arg_size == 8, "unhandled float size: %u", arg_size);
            src = Imm_Locs[dst_fn_arg_i >= 4][arg_size/4 - 1][is_float];
        }
        else
        {   src = Arg_Locs[src_fn_arg_i-- >= 4][is_float];   }

        // MOV AS APPROPRIATE //////////////////////////////////////////////////////////
        assert(dst && src);
        mc_size += Movs[dst][src](&cur, dsts[dst], srcs[src], Bnc_ax); // NOTE: ax is only used in mem_mem, where it's the tmp register
    }

    if (args_n > 4)
    { // NOTE: we have changed rsp, so we need control to return here to reset it
        mc_size += call_fn(&cur, fn);
        mc_size += mov_mem_mem(&cur, rsp_d, (BncVal){0}, Bnc_r11); // replace the return address // would it be better to just jmp to it?
        mc_size += add_rsp(&cur, rsp_d);
        mc_size += ret(&cur);
    }
    else
    {   mc_size += jmp_fn(&cur, fn);   }

    assert(!exe || (cur == (mem + data_size + mc_size)), "not calculating the pushed size correctly");
    return exe ?: (void *)(UPtr)(data_size + mc_size); // if we were passed a NULL, return the size, otherwise add nothing & return the start of the trampoline
}

#else // WIN32

#endif// WIN32


// mov
//   reg->reg
//     mov rax 4     b8 04 00 00 00
//     mov rcx 4     b9 04 00 00 00
//     mov rdx 4     ba 04 00 00 00
//     mov r8  4  41 b8 04 00 00 00
//     mov r9  4  41 b9 04 00 00 00
//   mem->reg
//   immediate->reg
//

static U8 args_3(U8 a, U8 b, U8 c)
{
    return a + b + c;
}

static U8 args_6(U8 a, U8 b, U8 c, U8 d, U8 e, U8 f)
{
    return a + b + c + d + e + f;
}

static U8 args_7(U8 a, U8 b, U8 c, U8 d, U8 e, U8 f, U8 g_)
{
    return a + b + c + d + e + f + g_;
}

static U8 args_16(U8 a, U8 b, U8 c, U8 d, U8 e, U8 f, U8 g_, U8 h,
                  U8 i, U8 j, U8 k, U8 l, U8 m, U8 n, U8  o, U8 p)
{
    return (a + b + c + d + e + f + g_ + h +
            i + j + k + l + m + n +  o + p);
}

int sub(int a, int b)
{ return a-b; }

int dec(int a)
{   return sub(a, 1);   }

typedef struct MyData {
    int Placeholder;
    float fs[6];
} MyData;

static float my_data_user(MyData data)
{
    return data.fs[3];
}

static Byte array_fn(Byte *arr, Size arr_n)
{
    return arr[arr_n-1];
}

static void spec_test()
{
    void *exe_buf = VirtualAlloc(NULL, 4096, MEM_COMMIT, PAGE_EXECUTE_READWRITE);
    U8 res_6 = args_6(1,2,3,4,0x123456789abcd,6);

    MyData data = { 42, .fs[3] = 6.28f };
    F4 (*fn2)() = specialize(exe_buf, my_data_user, ARRAY__N((BncArg[]) { bnc_struct(&data, sizeof(data)) }));
    F4 my_data_expect = my_data_user(data);
    F4 my_data_actual = fn2();

#define TEST(fmt, r, args, test_fn, array, expect, actual) do {\
    r (*fn) args = specialize(exe_buf, test_fn, ARRAY__N(array)); \
    r res[2] = { test_fn expect, fn actual }; \
    printf("%s "   fmt " = %s\n" \
           "     " fmt " = %s\n", \
           (res[0]==res[1] ? "PASS" : "FAIL"), \
           res[0], #test_fn #expect, \
           res[1], "fn"     #actual); \
    } while (0)

    Byte arr[3] = { 1, 2, 3 };
    Byte (*fn3)() = specialize(exe_buf, array_fn, ARRAY__N((BncArg[]) { bnc_copy(arr, sizeof(arr)), bnc_U8(sizeof(arr)) }));
    Byte byte = fn3();



    {
        U8 a=1, b=2, c=3, d=4, e=5, f=6, g_=7, h=8, i=9, j=10, k=11, l=12, m=13, n=14, o=15, p=16;

        TEST("%zu",U8,(U8), args_3,
             ((BncArg[]) { bnc_U8(a),  {sizeof(U8)}, bnc_U8(c) }),
             (a, b, c),
             (b));

        TEST("%zu",U8,(U8, U8, U8, U8), args_6,
             ((BncArg[]) { {sizeof(U8)}, {sizeof(U8)}, bnc_U8(c), {sizeof(U8)}, bnc_U8(e), {sizeof(U8)} }),
             (a, b, c, d, e, f),
             (a, b, d, f));

        TEST("%zu",U8,(U8, U8, U8, U8), args_7,
             ((BncArg[]) { {8}, {8}, bnc_U8(c), {8}, bnc_U8(e), {8}, bnc_U8(g_)}),
             (a, b, c, d, e, f, g_),
             (a, b, d, f));

        TEST("%zu",U8,(U8, U8, U8, U8), args_16,
             ((BncArg[]) {
              {Bnc_reg}, {Bnc_reg}, bnc_U8(c), {Bnc_reg}, bnc_U8(e), bnc_U8(f), bnc_U8(g_), bnc_U8(h),
              bnc_U8(i), bnc_U8(j), bnc_U8(k), bnc_U8(l), bnc_U8(m), {Bnc_reg}, bnc_U8(o),  bnc_U8(p), }),
             (a, b, c, d, e, f, g_, h,
              i, j, k, l, m, n, o,  p),
             (a, b, d, n));

        TEST("%zu",U8,(U8, U8, U8, U8), args_16,
             ((BncArg[]) {
              [2]=bnc_U8(c), [4]=bnc_U8(e), [5]=bnc_U8(f),  [6]=bnc_U8(g_), [7]=bnc_U8(h),
              [8]=bnc_U8(i), [9]=bnc_U8(j), [10]=bnc_U8(k), [11]=bnc_U8(l), [12]=bnc_U8(m), [14]=bnc_U8(o), [15]=bnc_U8(p), }),
             (a, b, c, d, e, f, g_, h,
              i, j, k, l, m, n, o,  p),
             (a, b, d, n));

        /* U8 (*fn2)(U8, U8, U8, U8) = specialize( */
        /*     exe_buf, args_16, ARRAY__N((BncArg[]) { */
        /*                                {BncArg_U8}, {BncArg_U8}, spec_U8(c), {BncArg_U8}, spec_U8(e), spec_U8(f),   spec_U8(g_), spec_U8(h), */
        /*                                spec_U8(i),   spec_U8(j),   spec_U8(k), spec_U8(l),   spec_U8(m), {BncArg_U8}, spec_U8(o),  spec_U8(p), })); */
        /* U8 res_16 = args_16(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16); */
        /* U8 res2 = fn2(a, b, d, n); */

        int stop = 0;
    }


    int stop2 = 0;
}
