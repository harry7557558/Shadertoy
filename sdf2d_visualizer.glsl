// 2D SDF Visualizer

// yellow-blue: positive/negative
// red-green: high/low gradient
// white: zero isoline

#define PI 3.1415926


float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) - k * h * (1.0 - h);
}
float smax(float a, float b, float k) {
    return -smin(-a, -b, k);
}
float sdf(in vec2 p) {
    if (p.y<-1.5) return p.y;
    if (p.y<-1.0) return -p.y;
    float r = 1.0 + 0.02*sin(16.0*atan(p.x-1.,p.y));
    float d1 = length(p-vec2(1,0)) - r;
    float b = 0.8 + 0.04*sin(10.0*p.x)*sin(10.0*p.y);
    float d2 = max(abs(p.x+1.0), abs(p.y)) - b;
    return smin(d1, d2, 0.5);
}


// visualization parameters
#define BOX_RADIUS 2.5
#define ISOLINE_WIDTH 0.1
#define SHOW_GRAD 1
#define HIGH_GRAD_HIGHLIGHT 1.0
#define LOW_GRAD_HIGHLIGHT 0.5


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    float e = (2.0*BOX_RADIUS)/min(iResolution.x,iResolution.y);
    vec2 p = (fragCoord-iMouse.xy)*e;

    float d = sdf(p);

#if SHOW_GRAD
    vec2 grad = vec2(
        sdf(p+vec2(e,0))-sdf(p-vec2(e,0)),
        sdf(p+vec2(0,e))-sdf(p-vec2(0,e))) / (2.0*e);
    float red = 1.0-exp(-HIGH_GRAD_HIGHLIGHT*max(length(grad)-1.0,0.));
    float green = 1.0-exp(-LOW_GRAD_HIGHLIGHT*max(1.0-length(grad),0.));
#else
    float red = 0.0;
    float green = 0.0;
#endif

    vec3 col = d>0.0 ? vec3(1.0,0.8,0.2) : vec3(0.3,0.6,1.0);  // sign
    col *= 0.7+0.3*cos(2.0*PI*d/ISOLINE_WIDTH);  // isolines
    col *= 1.0-0.8*exp(-0.5*abs(d)/ISOLINE_WIDTH);  // fade near zero
    col = mix(col, vec3(1.0), clamp(2.0-6.0*abs(d/ISOLINE_WIDTH),0.,1.));  // zero isoline
    col *= mix(vec3(1.0), vec3(2,0,0), red);  // high grad
    col *= mix(vec3(1.0), vec3(0,2,0), green);  // low grad
    fragColor = vec4(col, 1.0);
}
