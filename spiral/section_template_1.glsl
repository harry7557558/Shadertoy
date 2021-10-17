// cross section of a nautilus/spirula/snail shell

#define PI 3.1415926


const float r0 = 0.7;
const float h0 = 1.5;
const float b = 0.8;

// contains discontinuities
float sdf_min(in vec2 p) {
    // (x-h(i))**2 + y**2 == h(i)**2
    float h = dot(p,p)/(2.0*p.x);
    float i = log(h/h0)/b;
    float i0 = floor(i);
    //float min_d = 1e12;
    float min_d = length(p);
    for (float di=-0.0; di<=2.0; di++) {
        float i = i0+di;
        float h = h0*exp(b*i);
        float r = r0*exp(b*i);
        float d = length(p-vec2(h,0))-r;
        //min_d = min(min_d, d);
        if (min_d>0.0) min_d = d>0. ? min(d, min_d) : max(d, -min_d);
    }
    return min_d;
    return abs(min_d)-0.2;
}

float sdf(in vec2 p) {
    return sdf_min(p);
    return length(p)-1.0;
}


#define BOX_RADIUS 4.0
#define ISOLINE_WIDTH 0.1

#define SHOW_GRAD 1
#define HIGH_GRAD_HIGHLIGHT 1.0


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    float e = (2.0*BOX_RADIUS)/min(iResolution.x,iResolution.y);
    vec2 p = (fragCoord-iMouse.xy)*e;

    float d = sdf(p);

#if SHOW_GRAD
    vec2 grad = vec2(
        sdf(p+vec2(e,0))-sdf(p-vec2(e,0)),
        sdf(p+vec2(0,e))-sdf(p-vec2(0,e))) / (2.0*e);
    float red = 1.0-exp(-HIGH_GRAD_HIGHLIGHT*max(length(grad)-1.0,0.));
#else
    float red = 0.0;
#endif

    vec3 col = d>0.0 ? vec3(1.0,0.8,0.2) : vec3(0.3,0.6,1.0);
    col = mix(col, vec3(1,0,0), red);
    col *= 0.7+0.3*cos(2.0*PI*d/ISOLINE_WIDTH);
    col *= 1.0-0.8*exp(-0.5*abs(d)/ISOLINE_WIDTH);

    col = mix(col, vec3(1.0), clamp(2.0-6.0*abs(d/ISOLINE_WIDTH),0.,1.));

    fragColor = vec4(col, 1.0);
}
