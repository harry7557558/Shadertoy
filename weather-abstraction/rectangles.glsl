#include "noise.glsl"

#iChannel0 "img1.png"
#iChannel1 "img1_blurred.png"


void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    vec2 uv = fragCoord/iResolution.xy;
    vec2 b0 = vec2(0), b1 = vec2(1);

    float n = 8.0*texture(iChannel1, uv).x;

    for (float i=0.; i<8.; i++) {
        float h = hash12(mat2(cos(i),-sin(i),sin(i),cos(i))*0.5*(b0+b1))-0.5;
        h -= 0.5*(1.273*atan(b1.y-b0.y,b1.x-b0.x)-1.);
        if (h > 0.) {  // split into left and right
            float m = 0.5 + 0.2*(-1.0+2.0*hash12(b0+b1));
            float c = mix(b0.x, b1.x, m);
            if (uv.x < c) b1.x = c;
            else b0.x = c;
        }
        else {  // split into top and bottom
            float m = 0.5 + 0.2*(-1.0+2.0*hash12(b0+b1));
            float c = mix(b0.y, b1.y, m);
            if (uv.y < c) b1.y = c;
            else b0.y = c;
        }
        float v = length(texture(iChannel1, 0.5*(b0+b1)).xyz - vec3(.7,.6,.4));
        if (i > 3. && v > 0.5) break;
    }

    vec2 c = 0.5*(b0+b1);
    vec3 col = mix(texture(iChannel0,c), (
        texture(iChannel0,mix(c,vec2(b0.x,b0.y),0.5))+
        texture(iChannel0,mix(c,vec2(b0.x,b1.y),0.5))+
        texture(iChannel0,mix(c,vec2(b1.x,b0.y),0.5))+
        texture(iChannel0,mix(c,vec2(b1.x,b1.y),0.5)))/4., 0.4).xyz;
	col *= pow(16.0*c.x*c.y*(1.0-c.x)*(1.0-c.y), 0.2);
    col = 1.1*pow(col, vec3(1.4,1.4,1.4));

    fragColor = vec4(col, 1.0);
}
