#define PI 3.1415926

// http://iquilezles.org/www/articles/smin/smin.htm
float smin( float a, float b, float k )
{
    //return min(a, b);
    float h = max(k-abs(a-b),0.0);
    return min(a, b) - h*h*0.25/k;
}

// http://iquilezles.org/www/articles/smin/smin.htm
float smax( float a, float b, float k )
{
    //return max(a, b);
    float h = max(k-abs(a-b),0.0);
    return max(a, b) + h*h*0.25/k;
}

float sdEllipsoid( in vec3 p, in vec3 c, in vec3 r )
{
    return (length( (p-c)/r ) - 1.0) * min(min(r.x,r.y),r.z);
}
float sdTorus( vec3 p, vec2 t )
{
    return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}


float mapShell( in vec3 p ) 
{
    p -= vec3(0.05,0.12,-0.09);    

    vec3 q = mat3(-0.6333234236, -0.7332753384, 0.2474039592,
                   0.7738444477, -0.6034162289, 0.1924931824,
                   0.0081370606,  0.3133626215, 0.9495986813) * p;
    
    //q = p;

    const float b = 0.1759;
    
    float r = length( q.xy );
    float t = atan( q.y, q.x );
 
    // https://swiftcoder.wordpress.com/2010/06/21/logarithmic-spiral-distance-field/
    float np = (log(   r)/b-t)/(2.0*PI);
    float nm = (log(0.11)/b-t)/(2.0*PI);
    float n = min(np,nm);
    //n = np;
    
    float n0 = floor( n ), n1 = ceil(n);
    
    float r0 = exp( b * (t + 2.0*PI*n0));
    float r1 = exp( b * (t + 2.0*PI*n1));
    
    //-------

    float h0 = q.z + 1.5*r0 - 0.5;
    float h1 = q.z + 1.5*r1 - 0.5;
    //h0 = h1 = q.z;
    float d0 = sqrt((r0-r)*(r0-r)+h0*h0) - 0.99*r0;
    float d1 = sqrt((r1-r)*(r1-r)+h1*h1) - 0.99*r1;
    //d0 = 1e10;
    
    float d, dx, dy;
    if( d0<d1 ) { d = d0; dx=r0-r; dy=h0; }
    else        { d = d1; dx=r1-r; dy=h1; }
    //return 0.7*d;

    //float di = textureLod( iChannel2, vec2(t+r,0.5), 0. ).x;
    float di = 0.5*sin(40.*(t+r))+0.5*cos(60.*atan(dy,dx));
    d += 0.005*r*di;
    
    //matInfo = vec4(dx,dy,r/0.4,t/pi);

    vec3 s = q;
    q = q - vec3(0.34,-0.1,0.03);
    q.xy = mat2(0.8,0.6,-0.6,0.8)*q.xy;
    d = smin( d, sdTorus( q, vec2(0.28,0.05) ), 0.06);
    d = smax( d, -sdEllipsoid(q,vec3(0.0,0.0,0.0),vec3(0.24,0.36,0.24) ), 0.03 );
    d = smax( d, -sdEllipsoid(s,vec3(0.52,-0.0,0.0),vec3(0.42,0.23,0.5) ), 0.05 );
    
    return 0.7*d;
}


float map(in vec3 p) {
    return mapShell(p);
}

vec3 mapNormal(vec3 p) {
    mat3 k = mat3(p,p,p) - mat3(0.001);
    return normalize(map(p) - vec3(map(k[0]),map(k[1]),map(k[2])));
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    float rx = iMouse.z!=0.?3.14*(iMouse.y/iResolution.y)-1.57:0.3;
    float rz = iMouse.z!=0.?-iMouse.x/iResolution.x*4.0*3.14:0.5*iTime-2.0;

    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);

    vec3 ro = 3.0*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 2.0*length(iResolution.xy));
    rd = normalize(rd);

    float t0 = 0.01;
    float t1 = 2.0*length(ro);
    float t = t0;
    for (int i=min(int(iTime),0); i<100; i++) {
        float dt = map(ro+rd*t);
        t += dt;
        if (abs(dt) < 1e-3) break;
        if (t > t1) {
            fragColor = vec4(0, 0, 0, 1);
            return;
        }
    }
    vec3 light = normalize(vec3(0.5, 0.5, 1.0));
    //vec3 light = normalize(w+u+v);
    vec3 n = mapNormal(ro+rd*t);
    float col = 0.2+0.1*n.y+0.7*max(dot(n, light), 0.0);
    fragColor = vec4(vec3(col), 1.0);
    return;
}
