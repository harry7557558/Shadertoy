#include "noise.glsl"

#iChannel0 "img1.png"
#iChannel0::WrapMode "Clamp"


void test_triangle(vec2 p0, vec2 p1, vec2 p2, vec2 p, inout vec3 col) {
    if ((int(det(p0-p,p1-p)>0.)+int(det(p1-p,p2-p)>0.)+int(det(p2-p,p0-p)>0.))%3!=0) return;
    vec2 c = (p0+p1+p2)/3.;
    col = mix(texture(iChannel0,c), (
        texture(iChannel0,mix(c,p0,0.5))+
        texture(iChannel0,mix(c,p1,0.5))+
        texture(iChannel0,mix(c,p2,0.5)))/3., 0.4).xyz;
    //col = 0.1*floor(10.0*col);
    col = pow(col, vec3(1.4));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    vec2 sc = 0.18*iResolution.y/iResolution.xy;
    vec2 p = fragCoord.xy / iResolution.xy;
    vec2 pi = floor(p/sc+0.5)*sc;

    vec2 nbs[9];
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            vec2 q = pi + sc*(vec2(i,j)-1.);
            nbs[3*i+j] = q + 0.6*sc*(hash22(q)-0.5);
        }
    }

    vec3 col = vec3(0.0);
    for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
            vec2 t00 = nbs[3*i+j];
            vec2 t10 = nbs[3*(i+1)+j];
            vec2 t01 = nbs[3*i+(j+1)];
            vec2 t11 = nbs[3*(i+1)+(j+1)];
            if (length(t10-t01)<length(t11-t00)) {
                test_triangle(t00, t10, t01, p, col);
                test_triangle(t11, t01, t10, p, col);
            }
            else {
                test_triangle(t00, t11, t01, p, col);
                test_triangle(t11, t00, t10, p, col);
            }
        }
    }

    fragColor = vec4(col, 1.0);
}
