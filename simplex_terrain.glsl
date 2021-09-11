// hash functions from https://www.shadertoy.com/view/4djSRW by David Hoskins
float hash12(vec2 p) {
	vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}
float hash13(vec3 p3) {
	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.zyx + 31.32);
    return fract((p3.x + p3.y) * p3.z);
}
vec3 hash33(vec3 p3) {
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

// for solving interpolation
float det(vec2 a, vec2 b) {
    return a.x*b.y - a.y*b.x;
}
float det(vec3 a, vec3 b, vec3 c) {
    return dot(a, cross(b, c));
}

// noise functions and their analytical gradient

float SimplexValueNoise2D(vec2 xy) {
	// simplex grid
	const float K1 = 0.3660254038;  // (sqrt(3)-1)/2
	const float K2 = 0.2113248654;  // (-sqrt(3)+3)/6
	vec2 p = xy + (xy.x + xy.y)*K1;
	vec2 p1 = floor(p);
	vec2 s = xy.x-p1.x < xy.y-p1.y ? vec2(0., 1.) : vec2(1., 0.);
	vec2 p2 = p1 + s;
	vec2 p3 = p1 + 1.0;
	float v1 = 2.0 * hash12(p1) - 1.0;
	float v2 = 2.0 * hash12(p2) - 1.0;
	float v3 = 2.0 * hash12(p3) - 1.0;
	// interpolation
	vec2 f = p - p1, c = -s + 1.0;
	float m = 1.0 / det(s, c);
	float u = m * det(f, c);
	float uv = m * det(s, f);
	return v1 + u * (v2 - v1) + uv * (v3 - v2);  // mix(v1, mix(v2, v3, v), u)
}

vec3 SimplexValueNoise2Dg(vec2 xy) {
	// simplex grid
	const float K1 = 0.3660254038;  // (sqrt(3)-1)/2
	const float K2 = 0.2113248654;  // (-sqrt(3)+3)/6
	vec2 p = xy + (xy.x + xy.y)*K1;
	vec2 p1 = floor(p);
	vec2 s = xy.x-p1.x < xy.y-p1.y ? vec2(0., 1.) : vec2(1., 0.);
	vec2 p2 = p1 + s;
	vec2 p3 = p1 + 1.0;
	float v1 = 2.0 * hash12(p1) - 1.0;
	float v2 = 2.0 * hash12(p2) - 1.0;
	float v3 = 2.0 * hash12(p3) - 1.0;
	// interpolation
	vec2 f = p - p1, c = -s + 1.0;
	float m = 1.0 / det(s, c);
	float u = m * det(f, c);
	float uv = m * det(s, f);
	vec2 grad_u = m * vec2(c.y, -c.x);
	vec2 grad_uv = m * vec2(-s.y, s.x);
	float val = v1 + u * (v2 - v1) + uv * (v3 - v2);
	vec2 grad = grad_u * (v2 - v1) + grad_uv * (v3 - v2);
	return vec3(grad + (grad.x + grad.y)*K1, val);
}

float SimplexValueNoise3D(vec3 xyz) {
	const float K1 = 0.3333333333;
	const float K2 = 0.1666666667;
	// simplex grid
	vec3 p = xyz + (xyz.x + xyz.y + xyz.z)*K1;
	vec3 i = floor(p), f = p - i;
	vec3 f0 = xyz - (i - (i.x + i.y + i.z)*K2);
	//vec3 e = step(f0.yzx, f0);  // possibly result in degenerated simplex
	vec3 e = vec3(f0.y > f0.x ? 0.0 : 1.0, f0.z >= f0.y ? 0.0 : 1.0, f0.x > f0.z ? 0.0 : 1.0);
	vec3 i1 = e * (vec3(1.0) - e.zxy);
	vec3 i2 = vec3(1.0) - e.zxy * (vec3(1.0) - e);
	vec3 p0 = i;
	vec3 p1 = i + i1;
	vec3 p2 = i + i2;
	vec3 p3 = i + 1.0;
	float v0 = 2.0 * hash13(p0) - 1.0;
	float v1 = 2.0 * hash13(p1) - 1.0;
	float v2 = 2.0 * hash13(p2) - 1.0;
	float v3 = 2.0 * hash13(p3) - 1.0;
	// interpolation
	vec3 p01 = p1 - p0, p12 = p2 - p1, p23 = p3 - p2;
	float m = 1.0 / det(p01, p12, p23);
	float w = m * det(f, p12, p23);
	float uw = m * det(p01, f, p23);
	float uvw = m * det(p01, p12, f);
	return v0 + (v1 - v0) * w + (v2 - v1) * uw + (v3 - v2) * uvw;  // mix(v0, mix(mix(v1, v2, u), mix(v1, v3, u), v), w)
}

vec4 SimplexValueNoise3Dg(vec3 xyz) {
	const float K1 = 0.3333333333;
	const float K2 = 0.1666666667;
	// simplex grid
	vec3 p = xyz + (xyz.x + xyz.y + xyz.z)*K1;
	vec3 i = floor(p), f = p - i;
	vec3 f0 = xyz - (i - (i.x + i.y + i.z)*K2);
	vec3 e = vec3(f0.y > f0.x ? 0.0 : 1.0, f0.z >= f0.y ? 0.0 : 1.0, f0.x > f0.z ? 0.0 : 1.0);
	vec3 i1 = e * (vec3(1.0) - e.zxy);
	vec3 i2 = vec3(1.0) - e.zxy * (vec3(1.0) - e);
	vec3 p0 = i;
	vec3 p1 = i + i1;
	vec3 p2 = i + i2;
	vec3 p3 = i + 1.0;
	float v0 = 2.0 * hash13(p0) - 1.0;
	float v1 = 2.0 * hash13(p1) - 1.0;
	float v2 = 2.0 * hash13(p2) - 1.0;
	float v3 = 2.0 * hash13(p3) - 1.0;
	// interpolation
	vec3 p01 = p1 - p0, p12 = p2 - p1, p23 = p3 - p2;
	float m = 1.0 / det(p01, p12, p23);
	float w = m * det(f, p12, p23);
	float uw = m * det(p01, f, p23);
	float uvw = m * det(p01, p12, f);
	vec3 grad_w = m * cross(p12, p23);
	vec3 grad_uw = m * cross(p23, p01);
	vec3 grad_uvw = m * cross(p01, p12);
	float val = v0 + (v1 - v0) * w + (v2 - v1) * uw + (v3 - v2) * uvw;
	vec3 grad = (v1 - v0) * grad_w + (v2 - v1) * grad_uw + (v3 - v2) * grad_uvw;
	return vec4(grad + (grad.x + grad.y + grad.z)*K1, val);
}


// scene
float map(in vec3 p) {
    float base = p.z - 0.1*SimplexValueNoise2D(p.xy);
    float top = p.z - (0.05 * 2.0*SimplexValueNoise2D(0.5*p.xy) + 3.0);
    float blob = 0.6*5.0*SimplexValueNoise3D(0.2*p) + 0.2*SimplexValueNoise3D(p);
    return min(max(blob, top), base);
}

// to show that analytical gradient works
vec3 mapNormal(vec3 p) {
    vec3 t = SimplexValueNoise2Dg(p.xy);
    vec4 base = vec4(0, 0, 1, p.z) - 0.1*vec4(t.xy, 0.0, t.z);
    t = 2.0 * vec3(0.5,0.5,1)*SimplexValueNoise2Dg(0.5*p.xy);
    vec4 top = vec4(0, 0, 1, p.z) - (0.05 * vec4(t.xy, 0.0, t.z) + vec4(vec3(0), 3.0));
    vec4 blob = 0.6*5.0*vec4(vec3(0.2),1)*SimplexValueNoise3Dg(0.2*p) + 0.2*SimplexValueNoise3Dg(p);
    vec4 r = blob.w > top.w ? blob : top;
    r = r.w < base.w ? r : base;
    return normalize(r.xyz);
}



#define ZERO min(int(iTime),0)

const vec3 sundir = normalize(vec3(0.5, 0.5, 1.0));

bool raymarch(vec3 ro, inout vec3 rd, inout float t) {
    float t0 = 0.01, t1 = 120.0;
    for (int i=ZERO; i<120; i++) {
        float dt = map(ro+rd*t);
        t += dt;
        if (abs(dt) < 1e-2) break;
        if (t > t1) {
            return false;
        }
        rd = normalize(rd + vec3(0, 0, .001*dt));
    }
    return true;
}

float calcShadow(vec3 ro) {
    float t = 0.1;
    for (int i=ZERO; i<20; i++) {
        float dt = map(ro+sundir*t);
        if (dt < 0.) return 0.1;
        t += max(dt, 0.1);
    }
    return 1.0;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    float rx = iMouse.z>0. ? 1.8*(iMouse.y/iResolution.y)-0.2 : 0.28;
    float rz = iMouse.z>0. ? -iMouse.x/iResolution.x*4.0*3.14 : 0.5+0.1*iTime;
    float ry = iMouse.z>0. ? 0.0 : 0.1;

    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = vec3(-sin(rx)*vec2(cos(rz),sin(rz)), cos(rx));
    vec3 w = cross(u, v);
    u = cos(ry)*u + sin(ry)*v;
    v = cross(w, u);

    vec3 ro = vec3(2.0*iTime, 0, 2.5) + 14.0*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 2.0*length(iResolution.xy));
    rd = normalize(rd);

    float t = 0.0;
    vec3 col;
    if (raymarch(ro, rd, t)) {
        ro += rd * t;
        vec3 n = mapNormal(ro);
        vec3 basecol = 0.5+0.4*(1.0-abs(n.z))*hash33(n);
        float shadow = calcShadow(ro);
        vec3 sunlight = basecol * shadow * max(dot(n, sundir), 0.0) * vec3(1.2,1.1,0.8);
        vec3 skylight = basecol * max(n.z, 0.0) * vec3(0.6,0.6,0.8);
        vec3 indirect = basecol * max(-dot(n, sundir), 0.0) * vec3(0.6) + 0.2*max(-n.z, 0.0);
        col = sunlight + skylight + indirect;
    }
    else {
        t = 1e12;
    }
    col = mix(vec3(0.5, 0.6, 0.7)-0.3*max(rd.z, 0.0), col, exp(-0.015*t));
    col += vec3(0.8, 0.6, 0.4) * pow(max(dot(rd, sundir), 0.0), 1.5);
    col = 1.1*pow(col, vec3(1.2));
    fragColor = vec4(vec3(col), 1.0);
    return;
}
