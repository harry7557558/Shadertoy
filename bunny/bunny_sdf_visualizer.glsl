// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

#define PI 3.1415926


// a debugging SDF, heavily inspired by iq's Snail
float bunnySDF(in vec3 p) {
    //if (dot(p, p) > 1.1)  return length(p) - 1.0;
    vec4 f00=sin(p.x*vec4(1.16,.66,2.11,.4)+p.y*vec4(.8,-.43,.75,-.33)+p.z*vec4(.88,.16,1.82,-.31)+vec4(-.55,-.2,.35,.12));
    vec4 f01=sin(p.x*vec4(.24,-.36,.37,.48)+p.y*vec4(2.27,.12,2.46,.33)+p.z*vec4(.61,-.95,-1.49,-.56)+vec4(.53,.12,.44,.33));
    vec4 f02=sin(p.x*vec4(-1.84,-.42,-1.68,.8)+p.y*vec4(.04,.33,.98,1.13)+p.z*vec4(.53,.6,-1.23,.2)+vec4(-.55,.08,-.26,-.39));
    vec4 f03=sin(p.x*vec4(-.53,-.13,-1.51,1.13)+p.y*vec4(.25,.62,.8,-.34)+p.z*vec4(-.87,-.91,.53,-1.25)+vec4(-.37,.07,-.25,-.95));
    vec4 f10=sin(mat4(.66,-.74,-.47,.31,-.33,-.15,.27,.53,1.05,-.58,-.7,.48,-.05,.28,.22,-.16)*f00+mat4(.2,-1.01,-.89,-.26,.24,.05,-.06,.18,.05,-.63,-.96,-.16,.15,-.6,.31,.2)*f01+mat4(-.41,.02,.27,-.79,-.09,-.32,.29,-.36,-.52,.12,-.07,-.61,0,-.01,-.39,.56)*f02+mat4(-.03,.35,-.44,-.09,-.55,-.44,-.5,-.43,-.56,0,.02,-.42,-.06,.32,.06,.82)*f03+vec4(.19,-.17,.04,-.42));
    vec4 f11=sin(mat4(.51,-.24,-.31,-.21,.48,-.46,-.35,.43,.09,.22,-.43,-.04,-.24,-.31,.07,.01)*f00+mat4(.41,.65,-.87,-.16,-.63,-.29,.16,-.3,.16,-1.03,-.56,.5,-.37,-.45,.21,-.12)*f01+mat4(.39,.93,-.17,.05,.04,.78,-.19,.03,-.7,-1.16,-.18,-.3,.18,-.5,-.03,-.28)*f02+mat4(.17,-.43,.31,-.57,.1,-.1,-.24,-.27,-.35,.87,.26,-.61,.3,-.84,.08,-.03)*f03+vec4(-.47,-.09,.36,-.18));
    vec4 f12=sin(mat4(-.5,.43,-.34,-.73,-.09,.22,.21,.37,.68,-.2,-.58,-.32,.54,.47,-.09,.63)*f00+mat4(-.27,-.01,.52,-.48,-.19,.02,.7,.37,.45,-.53,.58,-.67,.33,-.22,-.2,.4)*f01+mat4(-.5,.44,-.25,.47,-.15,-.21,-.61,-.61,.17,-.2,-.32,.39,.43,.2,-.34,.04)*f02+mat4(-.06,.11,.21,-.34,.15,-.02,.58,-.25,-.39,-.3,-.2,-.41,.13,.61,-.3,-.23)*f03+vec4(-.04,-.07,-.06,.1));
    vec4 f13=sin(mat4(.17,.09,-.31,.28,-.17,.46,.27,-.16,.15,-.04,-.25,.12,-.07,.25,-.11,-.37)*f00+mat4(-.46,-.23,.2,-.62,.11,-.07,-.31,-.31,.21,.34,.14,-.86,-.11,-.49,-.43,-.55)*f01+mat4(.03,-.01,.1,-.14,.19,-.09,.11,-.25,-.25,.39,-.04,-.18,.11,.22,.06,-.02)*f02+mat4(.07,-.1,-.22,.14,.11,-.09,-.6,-.37,-.09,-.39,.23,-.95,.55,-.42,-.03,-.46)*f03+vec4(.22,-.15,.06,-.26));
    vec4 f20=sin(mat4(-.14,-.14,.18,.22,.12,-.21,.66,.12,-.26,.09,-.21,.2,.02,-.08,-.06,.17)*f10+mat4(.08,.13,-.75,.48,.35,-.5,.75,-.49,-.67,-.46,-.03,-.43,-.11,.1,0,-.29)*f11+mat4(.26,-.37,-.08,.27,-.63,-.43,.14,-.31,-.32,-.49,.08,.2,-.4,-.53,-.33,.07)*f12+mat4(-.31,-.06,-.46,-.49,-.17,.43,.28,.11,-.53,.1,.38,.65,.28,.05,.23,-.18)*f13+vec4(-.15,.08,-.14,-.22));
    vec4 f21=sin(mat4(-.06,-.01,.29,-.02,-.4,-.27,.06,.01,.64,.36,.84,-.66,-.01,-.59,.08,.54)*f10+mat4(.47,.35,.34,.34,-.54,.45,-.25,.43,-.27,-.06,-.49,-.06,.76,-.19,-.22,-.34)*f11+mat4(-.37,.21,-.65,.54,-.15,-.15,.02,.03,-.3,.18,-.4,-.52,.05,.22,-.22,-.56)*f12+mat4(.51,-.53,.35,.35,-.23,-.17,-.1,-.37,-.16,.16,.19,.02,-.11,-.08,-.22,.04)*f13+vec4(-.21,-.28,.04,.06));
    vec4 f22=sin(mat4(.2,-.17,-.08,.3,.12,-.09,-.16,.38,.34,.21,-.4,-.28,-.26,-.67,.51,.29)*f10+mat4(.25,-.51,.33,-.35,-.35,.81,-.53,.81,.1,.29,.07,-.15,.5,-.79,-.12,.04)*f11+mat4(-.65,.05,.19,.18,.33,-.41,-.06,.49,.12,.83,-.77,.28,-.16,.03,-.64,.74)*f12+mat4(-.37,.11,-.37,.07,-.21,-1.17,.34,.23,.06,-.48,.15,-.19,.74,.17,.17,-.15)*f13+vec4(.35,.76,-.03,.29));
    vec4 f23=sin(mat4(.05,-.33,-.26,-.31,.09,.01,.15,.24,-.4,-.03,-.37,-.4,.23,.11,.16,.05)*f10+mat4(-.3,-.37,.27,.43,.06,.81,-.68,-.22,-.01,-.6,.24,.43,-.31,.26,.33,-.08)*f11+mat4(.03,.2,-.01,.18,-.01,.18,-.04,.51,.24,-.45,.07,.37,-.15,.34,.27,-.14)*f12+mat4(.28,.08,-.11,-.14,.54,.28,-.22,.19,-.16,-.48,-.38,.04,-.57,-.58,.01,-.46)*f13+vec4(-.48,-.03,.01,0));
    return dot(vec4(-.41,.83,-.47,-.53),f20)+dot(vec4(-.5,-.46,.33,.32),f21)+dot(vec4(-.44,-.22,-.29,.5),f22)+dot(vec4(-.66,-.55,-.83,.62),f23)+.2;
}
float sdf(in vec3 p) {
    return bunnySDF(p)-0.01;
    //return bunnySDF(p)+0.3;
}

vec3 sdfGrad(in vec3 p, in float e) {
	float a = sdf(p+vec3(e,e,e));
	float b = sdf(p+vec3(e,-e,-e));
	float c = sdf(p+vec3(-e,e,-e));
	float d = sdf(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}



// raymarching parameters
const vec3 BoxRadius = vec3(1.2);
#define STEP 0.02
#define MAX_STEP 60.


// rendering parameters
#define FIELD_EMISSION 0.3
#define DISCONTINUITY_OPACITY 0.1
#define SURFACE_GRADIENT 10.0

// light direction as global variable
vec3 light = normalize(vec3(0.5,0.5,1.0));

// colormaps - https://www.shadertoy.com/view/NsSSRK
vec3 colorSdf(float t) {
  float r = .385+.619*t+.238*cos(4.903*t-2.61);
  float g = -5.491+.959*t+6.089*cos(.968*t-.329);
  float b = 1.107-.734*t+.172*cos(6.07*t-2.741);
  return clamp(vec3(r, g, b), 0.0, 1.0);
}
vec3 colorNormal(float t) {
  float r = .529-.054*t+.55*cos(5.498*t+2.779);
  float g = .21+.512*t+.622*cos(4.817*t-1.552);
  float b = .602-.212*t+.569*cos(5.266*t+2.861);
  return clamp(vec3(r, g, b), 0.0, 1.0);
}


// modified from a volume rendering demo
// https://github.com/harry7557558/Graphics/blob/master/raytracing/webgl_volume/fs-source.glsl
vec3 render(in vec3 ro, in vec3 rd, float t0, float t1) {
    float step_count = min(ceil((t1-t0)/STEP), MAX_STEP);
    float t = t0, dt = (t1-t0) / step_count;
    vec3 totcol = vec3(0.0);
    float totabs = 1.0;
    float v_old, v;
    for (t = t0; t < t1; t += dt) {
        v = sdf(ro+rd*t);
        vec3 col = colorSdf(0.5+0.5*sin(12.0*PI*v));
        //float grad = length(sdfGrad(ro+rd*t,dt));
        float grad = t==t0 ? 0.0 : abs(v-v_old)/dt;
        float grad_abs = (1.0-grad)/dt;
        col = mix(vec3(1,0,0), col, clamp(exp(grad_abs),0.0,1.0));
        float absorb = FIELD_EMISSION+DISCONTINUITY_OPACITY*max(-grad_abs,0.0);
        totabs *= exp(-absorb*dt);
        totcol += col*absorb*totabs*dt;
        if (v < 0.0) break;
        v_old = v;
    }
    if (v > 0.0) return totcol;
    for (int s = 0; s < 4; s += 1) {
        v_old = v;
        dt *= -0.5;
        for (int i = 0; i < 2; i++) {
            t += dt;
            v = sdf(ro+rd*t);
            if (v*v_old<0.0) break;
        }
    }
    vec3 grad = sdfGrad(ro+rd*t, 1e-3);
    vec3 col = colorNormal(0.5+0.5*tanh(SURFACE_GRADIENT*(0.5*length(grad)-0.5)));
    col = 0.2+0.05*grad.y+col*max(dot(normalize(grad), light),0.0);
    return totcol + col * totabs;
}


// ray intersection with a box
bool boxIntersection(vec3 ro, vec3 rd, out float tn, out float tf) {
    vec3 inv_rd = 1.0 / rd;
    vec3 n = inv_rd*(ro);
    vec3 k = abs(inv_rd)*BoxRadius;
    vec3 t1 = -n - k, t2 = -n + k;
    tn = max(max(t1.x, t1.y), t1.z);
    tf = min(min(t2.x, t2.y), t2.z);
    if (tn > tf) return false;
    return true;
}

// main
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // set camera
    float rx = iMouse.z>0. ? 3.14*(iMouse.y/iResolution.y)-1.57 : 0.3;
    float rz = iMouse.z>0. ? -iMouse.x/iResolution.x*4.0*3.14 : -0.6;

    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);

    vec3 ro = 8.0*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 2.0*length(iResolution.xy));
    rd = normalize(rd);

    // calculate pixel color
    light = normalize(w+0.5*u+0.1*v);

    float t0, t1;
    if (!boxIntersection(ro, rd, t0, t1)) {
        fragColor = vec4(vec3(0.0), 1.0);
        return;
    }
    vec3 col = render(ro, rd, t0, t1);;
    fragColor = vec4(col, 1.0);
}
