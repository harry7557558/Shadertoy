#include "noise.glsl"

#iChannel0 "img1.png"
#iChannel0::WrapMode "Clamp"


void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    vec2 sc = 0.18*iResolution.y/iResolution.xy;
    sc *= vec2(1.0, 0.8);
    vec2 p = fragCoord.xy / iResolution.xy;
    vec2 pi = floor(p/sc+0.5)*sc;

    vec2 nbs[9];
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            vec2 q = pi + sc*(vec2(i,j)-1.);
            vec2 uv = hash22(q);
            nbs[3*i+j] = q + 0.5*sc*vec2(cos(6.283*uv.x),cos(6.283*uv.x))*sqrt(uv.y);
        }
    }

    int min_id = -1;
    float min_dist = 100.0;
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            float d = length(p/sc-nbs[3*i+j]/sc);
            if (d<min_dist) min_dist=d, min_id=3*i+j;
        }
    }

    vec3 ctrcol = texture(iChannel0, nbs[min_id]).xyz;
    vec3 circol = vec3(0.);
    for (float i=0.; i<6.; i++) {
        vec2 dp = sc * vec2(cos(6.283*i/6.), sin(6.283*i/6.));
        circol += texture(iChannel0, nbs[min_id]+0.5*dp).xyz / 6.;
    }
    vec3 col = mix(ctrcol, circol, 0.2);
    col = 1.2*pow(col, vec3(1.4));

    fragColor = vec4(col, 1.0);
}
