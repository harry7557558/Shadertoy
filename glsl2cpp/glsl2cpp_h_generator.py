# Help creating `glsl2cpp.h`

def generate_vec_int_constructor(name: str, n: int):
    assert 2 <= n <= 4
    comps = ['x', 'y', 'z', 'w']
    lines = []
    for i in range(2**n):
        params = []
        inits = []
        for b in range(n):
            typename = "int" if (i >> b) & 1 == 1 else "float"
            compname = comps[b]
            params.append(f"const {typename} &{compname}")
            inits.append(f"{compname}((float){compname})")
        s = f"explicit {name}({', '.join(params)}) :{', '.join(inits)} {{}}"
        lines.append(s)
    return '\n'.join(lines)


def generate_swizzle(n: int):
    # vec2 yx() const { return vec2(y, x); }
    assert 2 <= n <= 4
    comps = ['x', 'y', 'z', 'w']
    lines = []
    for l in range(2, 5):
        for i in range(n**l):
            params = []
            for b in range(l):
                params.append(comps[(i//(n**b)) % n])
            s = f"vec{l} {''.join(params)}() const {{ return vec{l}({', '.join(params)}); }}"
            lines.append(s)
    return '\n'.join(lines)


def generate_univar_fun(name: str, n: int):
    # vec2 abs(const vec2 &a) { return vec2(abs(a.x), abs(a.y)); }
    assert 2 <= n <= 4
    comps = ['x', 'y', 'z', 'w']
    funs = ["saturate", "radians", "degrees",
            "sin", "cos", "tan", "asin", "acos", "atan",
            "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
            "exp", "log", "exp2", "log2", "sqrt", "inversesqrt",
            "abs", "sign", "floor", "ceil", "round", "trunc", "fract"]
    lines = []
    for fun in funs:
        params = []
        for i in range(n):
            params.append(f"{fun}(v.{comps[i]})")
        s = f"{name} {fun}(const {name} &v) {{ return {name}({', '.join(params)}); }}"
        lines.append(s)
    return '\n'.join(lines)


if __name__ == "__main__":
    print(generate_univar_fun("vec4", 4))
