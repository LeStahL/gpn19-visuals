#version 130

uniform float iTime;
uniform float iFFTWidth;
uniform float iScale;
uniform float iHighScale;
uniform float iNBeats;
uniform float iDial0;
uniform float iDial6;
uniform float iDial7;
uniform vec2 iResolution;
uniform sampler1D iFFT;
uniform float iNote;
uniform float iPressure;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void lfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}

void rot3(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

// Box sdf
void dbox(in vec2 x, in vec2 b, out float d)
{
    vec2 da = abs(x)-b;
    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}
/*
//distance to spline with parameter t
float dist2(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum dist3ance to spline
void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec2 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;
    vec3 ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

	//discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        ds = dist2(p0,p1,p2,x,ui.x+ui.y-tau);
        return;
    }
    
    //three dist3inct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    ds = min(
        dist2(p0,p1,p2,x, t.x),
        min(
            dist2(p0,p1,p2,x,t.y),
            dist2(p0,p1,p2,x,t.z)
        )
    );
}

void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}*/

void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)
{
    vec3 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
}

// Add sdfs
void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = mix(sda, sdb, step(sdb.x, sda.x));
}

vec2 ind;
void scene(in vec3 x, out vec2 sdf)
{
    x.y += mix(.3,5.,iDial0)*iTime;
    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));
    x.xy = R*x.xy;
    
    float d,
        size = mix(.5,.01,iDial6);
    vec2 x2 = mod(x.xy,size)-.5*size;
	
    ind = (x.xy - x2)/size;
    dbox(x2, .5*size*c.xx, d);
    zextrude(x.z, -d-.005, .05, d);
    d = max(x.z,d);
    d = abs(d);
    sdf = vec2(d,2.);
    
    float r;
    rand(ind-1.e2*floor(iTime), r);
    //lfnoise(12.*ind-1.*iTime, r);
    //r = .5+.5*r;
    if(r > .7)
    {
        dbox(x2, .5*size*c.xx, d);
        zextrude(x.z, -d-.02, (.1+.15*iScale)*(r-.7)/.3, d);
        stroke(d, .001, d);
        add(sdf, vec2(d,1.), sdf);
    }
}

void normal(in vec3 x, out vec3 n)
{
    const float dx = 1.e-4;
    vec2 s, na;
    
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void colorize(in vec2 x, out vec3 col)
{
    x.y += mix(.3,5.,iDial0)*iTime;
    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));
    x.xy = R*x.xy;
    
    float d,
        size = mix(.5,.01,iDial6),
        r;
    vec2 x2 = mod(x.xy,size)-.5*size;
    
    rand(ind-1.e2*floor(iTime), r);
    //lfnoise(12.*ind-1.*iTime, r);
    //r = .5+.5*r;
    col = mix(.14*c.xxx, .33*c.xxx, r);
    dbox(x2, .35*size*c.xx, d);
    if(r > .9)
    {
        col = mix(col, mix(c.xxy, c.xxx, .8), sm(d));
        stroke(d, .0025, d);
        col = mix(col, mix(c.xyy,c.xxx,.8), sm(d));
        stroke(d-.004, .002, d);
        col = mix(col, c.xyy, sm(d));
    }
	else if(r > .8)
    {
        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));
        stroke(d, .0025, d);
        col = mix(col, mix(.7*c.xxy,c.xxx,.8), sm(d));
        stroke(d-.004, .002, d);
        col = mix(col, .7*c.xxy, sm(d));
    }
    else if(r > .7)
    {
        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));
        stroke(d, .0025, d);
        col = mix(col, mix(mix(c.xxy, c.xyy, .5),c.xxx,.8), sm(d));
        stroke(d-.004, .002, d);
        col = mix(col, mix(c.xxy, c.xyy, .5), sm(d));
    }
    
    // Truchet
    /*
    float da = floor(4.*r)*pi/2.;
    R = mat2(cos(da), sin(da),-sin(da),cos(da));
    x2 = R * x2;
    if(r > .3)
    {
    	dspline2(x2,-.5*size*c.xy, c.yy, -.5*size*c.yx, d);
        dspline2(x2,.5*size*c.xy, c.yy, .5*size*c.yx, da);
        
    }
    else
    {
        dlinesegment(x2,-.5*size*c.xy, .5*size*c.xy, d);
        dlinesegment(x2,-.5*size*c.yx, .5*size*c.yx, da);
    }
    d = min(d, da);
    stroke(d, .001, d);
    col = mix(col,  .9*c.xyy, sm(d));
    stroke(d-.004, .002, d);
    col = mix(col, .0*c.xyy, sm(d));
*/
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), 
        s;
    vec3 col = c.yyy, 
        o = c.yzx,
        r = c.xyy, 
        u = normalize(c.yxx), 
        t = c.yyy, 
        dir,
        n,
        x;
    int N = 100,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.15)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        if(x.z<-.05)
        {
            col = .2*c.xxx;
            i = N;
            break;
        }
        d += min(s.x,5.e-3);
        //d += s.x;
    }
    
    if(i < N)
    {
        normal(x,n);
        
        if(s.y == 1.)
        {
            vec3 l = normalize(x+c.xzx);
            vec3 c1;
            
            float r;
		    rand(ind-1.e2*floor(iTime), r);
            if(r > .9)
                col = c.xyy;
            else if(r > .8)
                col = .7*c.xxy;
            else if(r > .7)
                col = mix(c.xxy, c.xyy, .5);
            float sc = clamp((r-.7)/.3,0.,1.);
            col = mix(mix(col, c.xxx, .1*sc), .4*c.xyy, sc);
            col = .3*col
                + .9*col * abs(dot(l,n))
                + 1.3*col * abs(pow(dot(reflect(-l,n),dir),3.));
            col = mix(col, c.xxx, .4);
            col *= col;
            
            d = -(o.z)/dir.z;
            x = o + d * dir;
            scene(x,s);
            l = normalize(x+c.xzx);
            colorize(x.xy, c1);
            n = c.yyx;
            
            c1 = .1*c1
                + .8*c1 * abs(dot(l,n))
                + c1 * abs(pow(dot(reflect(-l,n),dir),3.));
            col = mix(col, c1, .3);
        }
        else if(s.y == 2.)
        {
            vec3 l = normalize(x+c.xzx);
            float r;
            
            colorize(x.xy, col);
            col = .1*col
                + .8*col * abs(dot(l,n))
                + col * abs(pow(dot(reflect(-l,n),dir),3.));
        }
    }
    col += col;
    col *= col;
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
