// TODO: is there a way to do default args by this approach? (i.e. if the last arg is not provided it gets set to default)
// TODO: add a fixup arg type, where the instruction memory of the value be set is returned so that it can be altered by the caller
//          - this could default to 0 (don't need a larger input array) or could include passing a default value
//          - spec_copy could actually use this by only passing a size and having the caller fill out the data
//          - careful with thread safety etc
//          - or could overwrite all array values with their eventual location & size
// TODO: type checking
// TODO: manual tail-call optimization? https://ekins.space/2017/06/07/trampoline/

enum {
    ret = 0xc3,
};

typedef enum SpecXmm {
    Spec_xmm0,
    Spec_xmm1,
    Spec_xmm2,
    Spec_xmm3,
} SpecXmm;

typedef enum SpecReg {
    Spec_ax,     // rax, eax, ax
    Spec_cx,
    Spec_dx,
    Spec_r8 = 8, // r8,  r8d, r8w
    Spec_r9,
    Spec_r10,
    Spec_r11, // NOTE: volatile & unused for both Windows & Linux => good scratch register
} SpecReg;

// return an integer
//   exe mem
//   set rax
//   ret



/* typedef U8 FnA(); */
/* internal FnA * */
/* exe_ret_val(Byte *exe, U8 val) */
/* { */
/*     (void)val; */
/*     Byte instrs[] = { */
/*         /1* 0xb8, 0x04, 0x00, 0x00, 0x00, // mov rax 4 *1/ */
/*         mov_lit8_rcx(val), */
/*         mov_reg_reg(rax, rcx), */
/*         ret, */
/*     }; */
/*     memcpy(exe, instrs, countof(instrs)); */
/*     return (FnA *)exe; */
/* } */

internal U8 get_val() { return 0xcc00deadbeef00cc; }

/* typedef U8 FnA(); */
/* internal FnA * */
/* exe_jmp_fn(Byte *exe, U8 val) */
/* { */
/*     (void)val; */
/*     Byte instrs[] = { */
/*         mov_lit8_rax(val), */
/*         jmp_rax, */
/*     }; */
/*     memcpy(exe, instrs, countof(instrs)); */
/*     return (FnA *)exe; */
/* } */

internal U4 square(U4 val) { return val * val; }

/* typedef U4 FnA(); */
/* internal FnA * */
/* exe_param_1(Byte *exe, void *fn, U8 val) */
/* { */
/*     UPtr f = (UPtr)fn; */
/*     Byte instrs[] = { */
/*         mov_lit8_rcx(val), */
/*         mov_lit8_rax(f), */
/*         jmp_rax, */
/*     }; */
/*     memcpy(exe, instrs, countof(instrs)); */
/*     return (FnA *)exe; */
/* } */

typedef enum SpecArgTag {
    SpecArg_unknown = 0,

    SpecArg_ptr,
    SpecArg_U8,
    SpecArg_S8,
    SpecArg_F4,
    SpecArg_F8,

    SpecArg_copy,   // shouldn't be used directly
    SpecArg_struct, // shouldn't be used directly

    SpecArg_spec = 0x1000,

    SpecArg_spec_ptr,
    SpecArg_spec_U8,
    SpecArg_spec_S8,
    SpecArg_spec_F4,
    SpecArg_spec_F8,

    SpecArg_spec_copy,   // copies the data into the trampoline & passes the pointer to it
    SpecArg_spec_struct, // copies the data into the trampoline & passes the pointer to it if < 8 bytes

    SpecArg_Count,
} SpecArgTag;

typedef struct SpecArg {
    SpecArgTag tag;
    U4         size;
    union { void *ptr,*copy; U8 r_U8; S8 r_S8; F4 x_F4; F8 x_F8; };
} SpecArg;

global U1 const SpecArg_Sizes[] = {
    [SpecArg_ptr] = sizeof(void*),
    [SpecArg_U8]  = sizeof(U8),
    [SpecArg_S8]  = sizeof(S8),
    [SpecArg_F4]  = sizeof(F4),
    [SpecArg_F8]  = sizeof(F8),
};

internal inline SpecArg spec_arg(SpecArgTag tag) { return (SpecArg){ tag, SpecArg_Sizes[tag], }; }
internal inline SpecArg spec_copy(void   *v, Size size) { return (SpecArg){ SpecArg_spec_copy,   size, .copy = v }; }
internal inline SpecArg spec_struct(void *v, Size size) { return (SpecArg){ SpecArg_spec_struct, size, .copy = v }; }

internal inline SpecArg spec_ptr(void *v) { return (SpecArg){ SpecArg_spec_ptr, sizeof(void *), .ptr  = v }; }
internal inline SpecArg spec_U8(U8     v) { return (SpecArg){ SpecArg_spec_U8,  sizeof(U8),     .r_U8 = v }; }
internal inline SpecArg spec_S8(S8     v) { return (SpecArg){ SpecArg_spec_S8,  sizeof(S8),     .r_S8 = v }; }
internal inline SpecArg spec_F4(F4     v) { return (SpecArg){ SpecArg_spec_F4,  sizeof(F4),     .x_F4 = v }; }
internal inline SpecArg spec_F8(F8     v) { return (SpecArg){ SpecArg_spec_F8,  sizeof(F8),     .x_F8 = v }; }


#if 1  // MACHINE CODE
internal inline Size mc_cpy(Byte **cur, Byte const *mc, Size size)
{
    if (*cur)
    {   memcpy(*cur, mc, size), *cur += size;   }
    return size;
}

internal inline Size mov_reg_lit4(Byte **cur, SpecReg reg, U4 v)
{
    Size mc_size = ((reg >= Spec_r8) ? mc_cpy(cur, &(Byte){ 0x41 }, 1) : 0); // extended prefix

    Byte mc[] = {
        0xb8 + (reg % Spec_r8),                        // mov literal into r
        (v&255), (v>>8)&255, (v>>16)&255, (v>>24)&255, // LE 4-byte value
    };
    mc_size += mc_cpy(cur, mc, sizeof(mc));

    return mc_size;
}

internal inline Size mov_reg_lit8(Byte **cur, SpecReg reg, U8 v)
{
    Byte mc[] = {
        0x48 + (reg >= Spec_r8),                            // 8-byte modifier prefix
        0xb8 + (reg %  Spec_r8),                            // mov literal into r
        (v&255),     (v>>8)&255,  (v>>16)&255, (v>>24)&255, // LE 8-byte value
        (v>>32)&255, (v>>40)&255, (v>>48)&255, (v>>56)&255,
    };
    return mc_cpy(cur, mc, sizeof(mc));
}


internal inline Size mov_reg_reg(Byte **cur, SpecReg dst_reg, SpecReg src_reg)
{
    // NOTE: 4-byte version takes 2 bytes instead of 3 for e_x but not r8d,r9d
    // given that we can only swap these low registers max once per function, it's not really worth the complexity for 1 byte...
    Byte mc[] = {
        (0x48 + (src_reg >= Spec_r8) + 4*(dst_reg >= Spec_r8)), // 8-byte mod prefix
        0x89,                                                   // mov from register
        (0xc0 + 8*(src_reg%8) + (dst_reg%8)),                   // src & dst reg
    };
    return mc_cpy(cur, mc, sizeof(mc));
}

internal inline Size jmp_fn(Byte **cur, void *fn)
{
    Size mc_size = mov_reg_lit8(cur, Spec_ax, (U8)fn); // mov rax, fn
    if (cur)
    {   cur[0][0]=0xff, cur[0][1]=0xe0; *cur += 2;   } // jmp rax
    return mc_size + 2;
}

internal inline Size call_fn(Byte **cur, void *fn)
{
    Size mc_size = mov_reg_lit8(cur, Spec_ax, (U8)fn); // mov rax, fn
    if (cur)
    {   cur[0][0]=0xff, cur[0][1]=0xd0; *cur += 2;   } // call rax
    return mc_size + 2;
}

internal inline Size mov_xmm_xmm(Byte **cur, SpecXmm dst_xmm, SpecXmm src_xmm)
{
    Byte mc[] = {
        0xf3, 0x0f,                   // prefix
        0x7e,                         // mov to xmm
        (0xc0 + 8*dst_xmm + src_xmm), // src & dst
    };
    return mc_cpy(cur, mc, sizeof(mc));
}

// TODO: there may be a better approach; this is quite a large instruction
internal inline Size mov_xmm_F4(Byte **cur, SpecXmm xmm, F4 v)
{
    Size mc_size = mov_reg_lit4(cur, Spec_ax, (union{F4 f; U4 u;}){v}.u);
    Byte mc[]    = {
        0x66, 0x0f,       // prefix
        0x6e,             // mov to xmm
        (0xc0 + 8*(xmm)), // dst reg
    };
    return mc_size + mc_cpy(cur, mc, sizeof(mc));
}

// TODO: there may be a better approach; this is quite a large instruction
internal inline Size mov_xmm_F8(Byte **cur, SpecXmm xmm, F8 v)
{
    Size mc_size = mov_reg_lit8(cur, Spec_ax, (union{F8 f; U8 u;}){v}.u);
    Byte mc[]    = {
        0x66, 0x48, 0x0f, // prefix
        0x6e,             // mov to xmm
        (0xc0 + 8*(xmm)), // dst reg
    };
    return mc_size + mc_cpy(cur, mc, sizeof(mc));
}

typedef enum { reg_mem, mem_reg } SpecMemRegDir;
internal inline Size mov_reg_mem_reg(Byte **cur, SpecMemRegDir dir, SpecReg reg, Size rsp_offset)
{
    Byte dir_mc[] = { [reg_mem]=0x8b, [mem_reg]=0x89, };
    if (rsp_offset < 0x7f)
    {
        Byte mc[] = {
            0x48 + 4*(reg >= Spec_r8),                   // prefix
            dir_mc[dir], 0x44 + 8*(reg % Spec_r8), 0x24, // mov to reg
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
internal inline Size mov_reg_mem(Byte **cur, SpecReg reg, Size rsp_offset) { return mov_reg_mem_reg(cur, reg_mem, reg, rsp_offset); }
internal inline Size mov_mem_reg(Byte **cur, Size rsp_offset, SpecReg reg) { return mov_reg_mem_reg(cur, mem_reg, reg, rsp_offset); }
internal inline Size mov_mem_mem(Byte **cur, SpecReg tmp_reg, Size dst_rsp_offset, Size src_rsp_offset)
{
    Size mc_size  = mov_reg_mem_reg(cur, reg_mem, tmp_reg, src_rsp_offset);
    mc_size      += mov_reg_mem_reg(cur, mem_reg, tmp_reg, dst_rsp_offset);
    return mc_size;
}

internal inline Size mov_mem_xmm(Byte **cur, Size rsp_offset, SpecXmm xmm)
{
    if (rsp_offset < 0x7f)
    {
        Byte mc[] = {
            0x66, 0x0f,         // prefix
            0xd6, 0x44 + 8*xmm, // mov from xmm
            0x24,
            (S1)rsp_offset      // to [rsp + rsp_offset]
        }; // add rsp -size
        return mc_cpy(cur, mc, sizeof(mc));
    }
    else
    {
        TODO("handle larger rsp offsets (%zu)", rsp_offset);
        return 0;
    }
}


internal inline Size mov_mem_lit8(Byte **cur, Size rsp_offset, U8 v)
{
    // TODO (perf,mem): can do much more succinctly with values that fit in 4 bytes
    Size mc_size = mov_reg_lit8(cur, Spec_ax, v);
    mc_size += mov_mem_reg(cur, rsp_offset, Spec_ax);
    return mc_size;
}

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

#endif // MACHINE CODE

internal inline Bool arg_needs_copy(SpecArg const *arg)
{
    return (arg->tag == SpecArg_spec_copy ||
            (arg->tag == SpecArg_spec_struct &&
             arg->size > sizeof(void*)));
}

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
specialize(Byte *mem, void *fn, SpecArg const args[], Size args_n)
{
#define assume_args_are_sorted
    persist const SpecReg Arg_I_Regs[] = { Spec_cx, Spec_dx, Spec_r8, Spec_r9, };

    Byte *data      = mem;
    Size  mc_size   = 0,
          data_size = 0,
          data_used = 0;

    Size gen_fn_args_n = 0; // generated function
    for (Size i = 0; i < args_n; ++i) // count the number of generated args
    {
        SpecArg arg = args[i];
        gen_fn_args_n += !(arg.tag & SpecArg_spec);
        if (arg_needs_copy(&arg))
        {
            data_size = spec_align_up_offset(data, data_size, arg.size);
            data_size += args[i].size;
        }
    }
    Size gen_fn_arg_i = gen_fn_args_n - 1;

    Byte *exe = (mem ? mem + data_size : 0); // make space for data before the executable
    Byte *cur = exe;

    Size rsp_d = 0;
    if (args_n > 4)
    {
        rsp_d    = args_n - gen_fn_args_n;
        rsp_d   += rsp_d & 1; // rsp 16 byte align
        rsp_d   *= sizeof(void*);
        mc_size += sub_rsp(&cur, rsp_d);
        mc_size += mov_mem_mem(&cur, Spec_ax, 0, rsp_d); // move the return address to rsp
    }

    // NOTE: need to make sure we don't clobber args before we use them elsewhere
    for (Size tgt_fn_arg_i = args_n; tgt_fn_arg_i-- > 0;)
    { // move args from a literal/generated position to their target position
        SpecArg arg     = args[tgt_fn_arg_i];
        SpecReg tgt_reg = Arg_I_Regs[tgt_fn_arg_i];
        SpecXmm tgt_xmm = tgt_fn_arg_i;
        Size    tgt_mem = 8*(tgt_fn_arg_i);


        if (! (arg.tag & SpecArg_spec))
        { // relocate the gen arg to the appropriate index for the target arg
            if (tgt_fn_arg_i != gen_fn_arg_i) // otherwise the arg data is already in the correct register
            {
                if (tgt_fn_arg_i < 4) switch (arg.tag) // both in & out params are in registers
                {
                    case SpecArg_ptr:
                    case SpecArg_S8:
                    case SpecArg_U8: mc_size += mov_reg_reg(&cur, tgt_reg, Arg_I_Regs[gen_fn_arg_i]); break;

                    case SpecArg_F4:
                    case SpecArg_F8: mc_size += mov_xmm_xmm(&cur, tgt_xmm, gen_fn_arg_i); break;

                    default: unreachable("unhandled arg tag: %d", arg.tag);
                }

                else if (gen_fn_arg_i < 4) switch (arg.tag)
                {
                    case SpecArg_ptr:
                    case SpecArg_S8:
                    case SpecArg_U8: mc_size += mov_mem_reg(&cur, tgt_mem, Arg_I_Regs[gen_fn_arg_i]); break;

                    case SpecArg_F4:
                    case SpecArg_F8: mc_size += mov_mem_xmm(&cur, tgt_mem, gen_fn_arg_i); break;

                    default: unreachable("unhandled arg tag: %d", arg.tag);
                }

                else
                {   mc_size += mov_mem_mem(&cur, Spec_ax, tgt_mem, 8*(gen_fn_arg_i));   }
            }
            --gen_fn_arg_i;
        }

        else if (tgt_fn_arg_i < 4) switch (arg.tag &~ SpecArg_spec)
        { // handle args that we are fixing in place/won't be visible as parameters to the generated fn
            case SpecArg_struct:
            case SpecArg_copy: {
                if (arg_needs_copy(&arg)) // too big to fit in a register
                { // push some data and pass a pointer to it
                    if (cur)
                    {   memcpy(&data[data_used], arg.copy, arg.size);   }
                    data_used = spec_align_up_offset(data, data_used, arg.size);
                    arg.ptr   = &data[data_used];
                    data_used += arg.size;
                }
                else // fits in a register; put it in one
                {   memcpy(&arg.r_U8, arg.copy, arg.size);   } // copy the data into the literal U8 that we push as a register
            } // fallthrough

            case SpecArg_ptr:
            case SpecArg_U8:
            case SpecArg_S8: mc_size += mov_reg_lit8(&cur, tgt_reg, arg.r_U8); break;
            case SpecArg_F4: mc_size += mov_xmm_F4(&cur,   tgt_xmm, arg.x_F4); break;
            case SpecArg_F8: mc_size += mov_xmm_F8(&cur,   tgt_xmm, arg.x_F8); break;

            default: unreachable("unhandled arg tag: %d", arg.tag);
        }

        else
        {
            assert(tgt_fn_arg_i >= 4);
            mc_size += mov_mem_lit8(&cur, 8*tgt_fn_arg_i, arg.r_U8);
        }
    }

    if (args_n <= 4)
    {   mc_size += jmp_fn(&cur, fn);   }
    else // NOTE: we have changed rsp, so we need control to return here to reset it
    {
        mc_size += call_fn(&cur, fn);
        mc_size += mov_mem_mem(&cur, Spec_cx, rsp_d, 0); // replace the return address // would it be better to just jmp to it?
        mc_size += add_rsp(&cur, rsp_d);
        if (cur) { *cur = 0xc3; ++cur; } ++mc_size; // ret
    }

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
    F4 (*fn2)() = specialize(exe_buf, my_data_user, ARRAY__N((SpecArg[]) { spec_struct(&data, sizeof(data)) }));
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
    Byte (*fn3)() = specialize(exe_buf, array_fn, ARRAY__N((SpecArg[]) { spec_copy(arr, sizeof(arr)), spec_U8(sizeof(arr)) }));
    Byte byte = fn3();



    {
        U8 a=1, b=2, c=3, d=4, e=5, f=6, g_=7, h=8, i=9, j=10, k=11, l=12, m=13, n=14, o=15, p=16;

        TEST("%zu",U8,(U8), args_3,
             ((SpecArg[]) { spec_U8(a),  {SpecArg_U8}, spec_U8(c) }),
             (a, b, c),
             (b));

        TEST("%zu",U8,(U8, U8, U8, U8), args_6,
             ((SpecArg[]) { {SpecArg_U8}, {SpecArg_U8}, spec_U8(c), {SpecArg_U8}, spec_U8(e), {SpecArg_U8} }),
             (a, b, c, d, e, f),
             (a, b, d, f));

        TEST("%zu",U8,(U8, U8, U8, U8), args_7,
             ((SpecArg[]) { {SpecArg_U8}, {SpecArg_U8}, spec_U8(c), {SpecArg_U8}, spec_U8(e), {SpecArg_U8}, spec_U8(g_)}),
             (a, b, c, d, e, f, g_),
             (a, b, d, f));

        TEST("%zu",U8,(U8, U8, U8, U8), args_16,
             ((SpecArg[]) {
              {SpecArg_U8}, {SpecArg_U8}, spec_U8(c), {SpecArg_U8}, spec_U8(e), spec_U8(f),   spec_U8(g_), spec_U8(h),
              spec_U8(i),   spec_U8(j),   spec_U8(k), spec_U8(l),   spec_U8(m), {SpecArg_U8}, spec_U8(o),  spec_U8(p), }),
             (a, b, c, d, e, f, g_, h,
              i, j, k, l, m, n, o,  p),
             (a, b, d, n));

        /* U8 (*fn2)(U8, U8, U8, U8) = specialize( */
        /*     exe_buf, args_16, ARRAY__N((SpecArg[]) { */
        /*                                {SpecArg_U8}, {SpecArg_U8}, spec_U8(c), {SpecArg_U8}, spec_U8(e), spec_U8(f),   spec_U8(g_), spec_U8(h), */
        /*                                spec_U8(i),   spec_U8(j),   spec_U8(k), spec_U8(l),   spec_U8(m), {SpecArg_U8}, spec_U8(o),  spec_U8(p), })); */
        /* U8 res_16 = args_16(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16); */
        /* U8 res2 = fn2(a, b, d, n); */

        int stop = 0;
    }


    int stop2 = 0;
}
