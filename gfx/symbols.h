//Generated with Symbolize (c) 2019 Alexander Kraus <nr4@z10.info>.
#ifndef SYMBOLIZE_H
#define SYMBOLIZE_H

int rand_handle, zextrude_handle, stroke_handle, smoothmin_handle, dhexagonpattern_handle, normal_handle, rot3_handle, lfnoise_handle, dbox_handle, dbox3_handle, dvoronoi_handle, dquadvoronoi_handle, analytical_box_handle, mfnoise_handle, dpolygon_handle, dstar_handle, dcirclesegment_handle, dcircle_handle, dlinesegment_handle, dlogo210_handle, rand3_handle, dspline3_handle;
const int nsymbols = 22;
const char *rand_source = "#version 130\n\n"
"void rand(in vec2 x, out float n)\n"
"{\n"
"    x += 400.;\n"
"    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\0";
const char *zextrude_source = "// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\0";
const char *stroke_source = "// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\0";
const char *smoothmin_source = "// iq's smooth minimum\n"
"void smoothmin(in float a, in float b, in float k, out float dst)\n"
"{\n"
"    float h = max( k-abs(a-b), 0.0 )/k;\n"
"    dst = min( a, b ) - h*h*h*k*(1.0/6.0);\n"
"}\n"
"\0";
const char *dhexagonpattern_source = "// Distance to hexagon pattern\n"
"void dhexagonpattern(in vec2 p, out float d, out vec2 ind) \n"
"{\n"
"    vec2 q = vec2( p.x*1.2, p.y + p.x*0.6 );\n"
"    \n"
"    vec2 pi = floor(q);\n"
"    vec2 pf = fract(q);\n"
"\n"
"    float v = mod(pi.x + pi.y, 3.0);\n"
"\n"
"    float ca = step(1.,v);\n"
"    float cb = step(2.,v);\n"
"    vec2  ma = step(pf.xy,pf.yx);\n"
"    \n"
"    d = dot( ma, 1.0-pf.yx + ca*(pf.x+pf.y-1.0) + cb*(pf.yx-2.0*pf.xy) );\n"
"    ind = pi + ca - cb*ma;\n"
"    ind = vec2(ind.x/1.2, ind.y);\n"
"    ind = vec2(ind.x, ind.y-ind.x*.6);\n"
"}\n"
"\0";
const char *normal_source = "const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"void scene(in vec3 x, out vec2 s);\n"
"void normal(in vec3 x, out vec3 n)\n"
"{\n"
"    const float dx = 5.e-4;\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\0";
const char *rot3_source = "const vec3 c = vec3(1.,0.,-1.);\n"
"void rot3(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\0";
const char *lfnoise_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void rand(in vec2 x, out float d);\n"
"void lfnoise(in vec2 t, out float n)\n"
"{\n"
"    vec2 i = floor(t);\n"
"    t = fract(t);\n"
"    t = smoothstep(c.yy, c.xx, t);\n"
"    vec2 v1, v2;\n"
"    rand(i, v1.x);\n"
"    rand(i+c.xy, v1.y);\n"
"    rand(i+c.yx, v2.x);\n"
"    rand(i+c.xx, v2.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    n = mix(v1.x, v1.y, t.x);\n"
"}\n"
"\0";
const char *dbox_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"    vec2 da = abs(x)-b;\n"
"    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\0";
const char *dbox3_source = "#version 130\n\n"
"void dbox3(in vec3 x, in vec3 b, out float d)\n"
"{\n"
"  vec3 da = abs(x) - b;\n"
"  d = length(max(da,0.0))\n"
"         + min(max(da.x,max(da.y,da.z)),0.0);\n"
"}\n"
"\0";
const char *dvoronoi_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void rand(in vec2 x, out float d);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z)\n"
"{\n"
"    vec2 y = floor(x);\n"
"       float ret = 1.;\n"
"    vec2 pf=c.yy, p;\n"
"    float df=10.;\n"
"    \n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            d = length(x-p);\n"
"            \n"
"            if(d < df)\n"
"            {\n"
"                df = d;\n"
"                pf = p;\n"
"            }\n"
"        }\n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            vec2 o = p - pf;\n"
"            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);\n"
"            ret = min(ret, d);\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    z = pf;\n"
"}\n"
"\0";
const char *dquadvoronoi_source = "#version 130\n\n"
"// Random Quadtree\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void rand(in vec2 x, out float num);\n"
"void dquadvoronoi(in vec2 x, in float threshold, in float depth, out float d, out float faco)\n"
"{\n"
"    d = 1.;\n"
"    vec2 y = x, \n"
"        yi;\n"
"    float size = .5,\n"
"	    fac = 1.;\n"
"    faco = 1.;\n"
"    for(float i=0.; i<depth; i+=1.)\n"
"    {\n"
"        float dd;\n"
"            vec2 ind;\n"
"            dvoronoi(y/size/.5,dd, ind);\n"
"        vec2 y0 = y;\n"
"		float r;\n"
"        rand(ind+fac,r);\n"
"        fac *= r*step(r,threshold);\n"
"        faco *= r;\n"
"        if(fac != 0.)\n"
"        {\n"
"            \n"
"            //dd = mix(dd,length(y)-.5*size,step(r,threshold));\n"
"            dd = abs(dd);\n"
"            smoothmin(d,dd,.01,d);\n"
"        }\n"
"        \n"
"        size *= .5;\n"
"    }\n"
"    faco += fac*fac;\n"
"}\n"
"\0";
const char *analytical_box_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void analytical_box(in vec3 o, in vec3 dir, in vec3 size, out float d)\n"
"{\n"
"    vec3 tlo = min((size-o)/dir,(-size-o)/dir); // Select 3 visible planes\n"
"    vec2 abxlo = abs(o.yz + tlo.x*dir.yz),\n"
"        abylo = abs(o.xz + tlo.y*dir.xz),\n"
"        abzlo = abs(o.xy + tlo.z*dir.xy);\n"
"    vec4 dn = 100.*c.xyyy;\n"
"    \n"
"    dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));\n"
"    dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));\n"
"    dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));\n"
"\n"
"    d = dn.r;\n"
"}\n"
"\0";
const char *mfnoise_source = "#version 130\n\n"
"// const vec3 c = vec3(1.,0.,-1.);\n"
"void lfnoise(in vec2 x, out float d);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)\n"
"{\n"
"    n = 0.;\n"
"    float a = 1., nf = 0., buf;\n"
"    for(float f = d; f<b; f *= 2.)\n"
"    {\n"
"        lfnoise(f*x, buf);\n"
"        n += a*buf;\n"
"        a *= e;\n"
"        nf += 1.;\n"
"    }\n"
"    n *= (1.-e)/(1.-pow(e, nf));\n"
"}\n"
"\0";
const char *dpolygon_source = "#version 130\n\n"
"// compute distance to regular polygon\n"
"const float pi = acos(-1.);\n"
"void dpolygon(in vec2 x, in float N, in float R, out float dst)\n"
"{\n"
"    float d = 2.*pi/N,\n"
"        t = mod(acos(x.x/length(x)), d)-.5*d;\n"
"    dst = R-length(x)*cos(t)/cos(.5*d);\n"
"}\n"
"\0";
const char *dstar_source = "#version 130\n\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void dstar(in vec2 x, in float N, in vec2 R, out float dst)\n"
"{\n"
"    float d = pi/N,\n"
"        p0 = acos(x.x/length(x)),\n"
"        p = mod(p0, d),\n"
"        i = mod(round((p-p0)/d),2.);\n"
"    x = length(x)*vec2(cos(p),sin(p));\n"
"    vec2 a = mix(R,R.yx,i),\n"
"    	p1 = a.x*c.xy,\n"
"        ff = a.y*vec2(cos(d),sin(d))-p1;\n"
"   	ff = ff.yx*c.zx;\n"
"    dst = dot(x-p1,ff)/length(ff);\n"
"}\n"
"\0";
const char *dcirclesegment_source = "#version 130\n\n"
"const float pi = acos(-1.);\n"
"void dcirclesegment(in vec2 x, in float R, in float p0, in float p1, out float d)\n"
"{\n"
"    float p = atan(x.y, x.x);\n"
"    vec2 philo = vec2(max(p0, p1), min(p0, p1));\n"
"    if((p < philo.x && p > philo.y) || (p+2.0*pi < philo.x && p+2.0*pi > philo.y) || (p-2.0*pi < philo.x && p-2.0*pi > philo.y))\n"
"        d = abs(length(x)-R);\n"
"    else d = min(\n"
"        length(x-vec2(cos(p0), sin(p0))),\n"
"        length(x-vec2(cos(p1), sin(p1)))\n"
"        );\n"
"}\n"
"\0";
const char *dcircle_source = "#version 130\n\n"
"void dcircle(in vec2 x, in float R, out float d)\n"
"{\n"
"    d = abs(length(x)-R);\n"
"}\n"
"\0";
const char *dlinesegment_source = "#version 130\n\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\0";
const char *dlogo210_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"void dcirclesegment(in vec2 x, in float R, in float p0, in float p1, out float d);\n"
"void dcircle(in vec2 x, in float R, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void dlogo210(in vec2 x, in float R, out float d)\n"
"{\n"
"    float d2;\n"
"    dcircle(x+R*c.zy, R, d);\n"
"    dlinesegment(x, c.yz, c.yx, d2);\n"
"    d = min(d, d2);\n"
"    dcirclesegment(x+R*c.xy, R, -.5*pi, .5*pi, d2);\n"
"    d = min(d, d2);\n"
"}\n"
"\0";
const char *rand3_source = "#version 130\n\n"
"void rand3(in vec3 x, out float num)\n"
"{\n"
"    x += 400.;\n"
"    num = fract(sin(dot(sign(x)*abs(x) ,vec3(12.9898,78.233,121.112)))*43758.5453);\n"
"}\n"
"\0";
const char *dspline3_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"\n"
"//distance to spline with parameter t\n"
"float dist(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)\n"
"{\n"
"    t = clamp(t, 0., 1.);\n"
"    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);\n"
"}\n"
"\n"
"//minimum distance to spline\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds)\n"
"{\n"
"    //coefficients for 0 = t^3 + a * t^2 + b * t + c\n"
"    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0,\n"
"    	ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);\n"
"\n"
"	//discriminant and helpers\n"
"    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;\n"
"    \n"
"    //triple real root\n"
"    if(dis > 0.) \n"
"    {\n"
"        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);\n"
"        ds = dist(p0,p1,p2,x,ui.x+ui.y-tau);\n"
"        return;\n"
"    }\n"
"    \n"
"    //three distinct real roots\n"
"    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;\n"
"    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;\n"
"    ds = min(\n"
"        dist(p0,p1,p2,x, t.x),\n"
"        min(\n"
"            dist(p0,p1,p2,x,t.y),\n"
"            dist(p0,p1,p2,x,t.z)\n"
"        )\n"
"    );\n"
"}\n"
"\0";
const char *hexagontunnel_source = "/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
"* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dhexagonpattern(in vec2 p, out float d, out vec2 ind);\n"
"void normal(in vec3 x, out vec3 n);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"void lfnoise(in vec2 t, out float n);\n"
"\n"
"float mat;\n"
"void scene(in vec3 x, out vec2 d)\n"
"{\n"
"    d = c.xx;\n"
"    \n"
"    x.z -= mix(1.,5.,iDial0)*iTime;\n"
"    \n"
"    float phi = atan(x.y, x.x),\n"
"        dhex,\n"
"        na,\n"
"        nal;\n"
"    vec2 ind;\n"
"    rand(floor(.33*iTime)*c.xx, na);\n"
"    rand(floor(.33*iTime)*c.xx+1., nal);\n"
"    na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));\n"
"    dhexagonpattern(mix(1.,4.,na)*1.01*vec2(pi,3.)*vec2(phi,x.z),dhex,ind);\n"
"    rand(ind,mat);\n"
"    stroke(dhex, .1, dhex);\n"
"    mat *= (1.-iScale);\n"
"    //d.x = length(x.xy)-mix(.5,.45,smoothstep(.1,-.1,dhex));\n"
"    d.x = min(d.x, length(x.xy)+.2*mat-mix(.5,.55+.2*mat,smoothstep(.2,-.2,dhex)));\n"
"    \n"
"    stroke(d.x, .04, d.x);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    // Set up coordinates\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(col, 0.);\n"
"        return;\n"
"    }\n"
"    \n"
"    // Camera setup\n"
"    float pp = .3*iTime;\n"
"    vec3 o = c.yyx,\n"
"        t = c.yyy,\n"
"        dir = normalize(t-o),\n"
"        r = normalize(c.xyy),\n"
"        u = normalize(cross(r,dir)),\n"
"        n,\n"
"        x,\n"
"        l;\n"
"    t += uv.x*r + uv.y*u;\n"
"    dir = normalize(t-o);\n"
"    vec2 s;\n"
"    float d = (.3+.2*mat)/length(dir.xy);// -(o.z-.12)/dir.z;\n"
"    int N = 450,\n"
"        i;\n"
"    \n"
"    // Graph\n"
"    x = o + d * dir;\n"
"    \n"
"    // Actual Scene\n"
"    {\n"
"\n"
"        // Raymarching\n"
"        for(i=0; i<N; ++i)\n"
"        {\n"
"            x = o + d * dir;\n"
"            scene(x,s);\n"
"            if(s.x < 1.e-4) break;\n"
"            d += min(s.x,.005);\n"
"        }\n"
"\n"
"        // Illumination\n"
"        l = normalize(x+c.yxx);\n"
"        if(i<N)\n"
"        {\n"
"            normal(x,n);\n"
"                        \n"
"            float phi = atan(x.y, x.x),\n"
"                dhex,\n"
"                na,\n"
"                nal;\n"
"            vec2 ind;\n"
"            rand(floor(.33*iTime)*c.xx, na);\n"
"            rand(floor(.33*iTime)*c.xx+1., nal);\n"
"            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));\n"
"            \n"
"            mat3 RR;\n"
"            rot3(na*1.e3*vec3(1.1,1.5,1.9)+12.*iNote,RR);\n"
"\n"
"            col = mix((.5+.5*mat)*c.xxx,(1.+.8*mat)*abs(RR*vec3(0.89,0.44,0.23)),.5+.5*sin(x.z));\n"
"            col = mix(col,vec3(0.25,0.23,0.21),.5+.5*cos(4.*x.z+mat));\n"
"\n"
"            \n"
"            dhexagonpattern(mix(1.,4.,na)*1.01*vec2(pi,3.)*vec2(phi,x.z-iTime),dhex,ind);\n"
"            stroke(dhex, .3, dhex);\n"
"            col = mix(col, clamp(1.9*col,c.yyy,c.xxx), mat*smoothstep(1.5/iResolution.y, -1.5/iResolution.y, -dhex));\n"
"        }\n"
"    }\n"
"    \n"
"    // Colorize\n"
"    col = .2*col\n"
"        + .9*col*abs(dot(l,n))\n"
"        +.7*col*abs(pow(dot(reflect(-l,n),dir),3.));\n"
"    \n"
"    fragColor = clamp(vec4(col,1.0),0.,1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *voronoinet_source = "/* Corfield Imitation 1\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" * \n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" * \n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" * \n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
" */\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dbox(in vec2 p, in vec2 b, out float dst);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void dquadvoronoi(in vec2 x, in float threshold, in float depth, out float d, out float faco);\n"
"void analytical_box(in vec3 o, in vec3 dir, in vec3 size, out float d);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"\n"
"// Scene\n"
"float mat;\n"
"void scene(in vec3 x, out vec2 d)\n"
"{\n"
"    d = c.xx;\n"
"    float dbound;\n"
"    dbox3(x,vec3(.3*c.xx, .2),dbound);\n"
"    float da, fac;\n"
"    dquadvoronoi(x.xy-.1*iTime, .71, 4., da, fac);\n"
"    \n"
"    float p = pi/4.;\n"
"    vec2 cs = vec2(cos(p),sin(p));\n"
"    mat2 m = mat2(cs.x,cs.y,-cs.y,cs.x);\n"
"    vec2 y = m*x.xy;\n"
"    float da9, fac9;\n"
"    dquadvoronoi(y-12.-.1*iTime, .41,2., da9, fac9);\n"
"    smoothmin(da,da9,.01,da);\n"
"    \n"
"    float r;\n"
"    rand(202.*fac*fac9*c.xx+3., r);\n"
"    mat = r;\n"
"    zextrude(x.z,da,r*.3,da9);\n"
"    smoothmin(d.x,da9,.4, d.x);\n"
"    \n"
"    stroke(da, .015+6.*.045*clamp(iScale,0.,1.), da);\n"
"    float db;\n"
"    stroke(da, .011+6.*.033*clamp(iScale,0.,1.), db);\n"
"   \n"
"    stroke(d.x,.003,d.x);\n"
"    dbox3(x,vec3(.33*c.xx, .02),da);\n"
"    smoothmin(d.x,da,.2,d.x);\n"
"    smoothmin(d.x,db,.05,d.x);\n"
"    //d.x = min(d.x, da);\n"
"    //d.x = min(d.x, db);\n"
"    \n"
"//     d.x = max(d.x, dbound);\n"
"}\n"
"\n"
"// Normal\n"
"const float dx = 5.e-4;\n"
"void normal(in vec3 x, out vec3 n);\n"
"\n"
"// Texture\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{    \n"
"    float phi = .1*iTime;\n"
"    \n"
"    vec3 white = .4*vec3(0.99,0.29,0.09),\n"
"        gray = vec3(0.95,0.25,0.05);\n"
"    float size = .1;\n"
"    \n"
"    \n"
"    vec2 y = mod(x,size)-.5*size;\n"
"    y = abs(y)-.001;\n"
"    float d = min(y.x,y.y);\n"
"\n"
"    y = mod(x,.2*size)-.1*size;\n"
"    y = abs(y)-.0002;\n"
"    d = min(d, min(y.x,y.y));\n"
"    \n"
"    col = mix(white, gray, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));\n"
"	col = mix(col, length(col)/length(c.xxx)*c.xxx, .7);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"     // Set up coordinates\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(col, 0.);\n"
"        return;\n"
"    }\n"
"    \n"
"    float dhex,\n"
"                na,\n"
"                nal;\n"
"            vec2 ind;\n"
"            rand(floor(.33*iTime)*c.xx, na);\n"
"            rand(floor(.33*iTime)*c.xx+1., nal);\n"
"            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));\n"
"    uv = mix(uv,.5*uv.yx,na);\n"
"    \n"
"    // Camera setup\n"
"    float pp = .3*iTime;\n"
"    vec3 o = c.yyx+.2*c.yzy,\n"
"        t = c.yyy,\n"
"        dir = normalize(t-o),\n"
"        r = normalize(c.xyy),\n"
"        u = normalize(cross(r,dir)),\n"
"        n,\n"
"        x,\n"
"        l;\n"
"    t += uv.x*r + uv.y*u;\n"
"    dir = normalize(t-o);\n"
"    vec2 s;\n"
"    float d = -(o.z-.2)/dir.z;\n"
"    int N = 250,\n"
"        i;\n"
"    \n"
"    // Graph\n"
"    //analytical_box(o,dir,vec3(.3*c.xx,.2),d);\n"
"    x = o + d * dir;\n"
"    \n"
"    // Actual Scene\n"
"    if(x.z>0.)\n"
"    {\n"
"\n"
"        // Raymarching\n"
"        for(i=0; i<N; ++i)\n"
"        {\n"
"            x = o + d * dir;\n"
"            scene(x,s);\n"
"            if(s.x < 1.e-4) break;\n"
"            d += min(.008,s.x);\n"
"        }\n"
"\n"
"        // Illumination\n"
"        l = normalize(x+c.yxx);\n"
"        if(i<N)\n"
"        {\n"
"            normal(x,n);\n"
"            \n"
"            mat3 RR;\n"
"            rot3(na*1.e3*vec3(1.1,1.5,1.9)+iNote*210.,RR);\n"
"            col = mix(mix(.0,.3,clamp(x.z/.3,0.,1.))*(.5+.5*mat)*c.xxx,(1.+.8*mat)*abs(RR*RR*vec3(.7,.5,.26)),step(x.z,.08));\n"
"            col = mix(col,(1.+.8*mat)*abs(RR*vec3(.6,.12,.06)),step(.19,x.z));\n"
"\n"
"            col = mix((.5+.5*mat)*col,(1.+.8*mat)*abs(RR*vec3(0.89,0.44,0.23)),(.5+.5*sin(x.z))*step(.19,x.z));\n"
"            col = mix(col,vec3(0.25,0.23,0.21),(.5+.5*cos(4.*x.z+mat))*step(.19,x.z));\n"
"            \n"
"            col = mix(col, clamp(1.9*col,c.yyy,c.xxx), mat*step(.19,x.z));\n"
"        }\n"
"        else\n"
"        {\n"
"            d = -o.z/dir.z;\n"
"            x = o + d * dir;\n"
"            n = c.yyx;\n"
"            l = vec3(x.xy, .8);\n"
"            colorize(x.xy, col);\n"
"        }\n"
"    }\n"
"    else // Floor with grid\n"
"    {\n"
"        d = -o.z/dir.z;\n"
"        x = o + d * dir;\n"
"        n = c.yyx;\n"
"        l = vec3(x.xy, .8);\n"
"        colorize(x.xy, col);\n"
"    }\n"
"    \n"
"    // Colorize\n"
"    col = .2*col\n"
"        + .9*col*abs(dot(l,n))\n"
"        +.4*col*abs(pow(dot(reflect(-l,n),dir),3.));\n"
"    \n"
"    fragColor = clamp(vec4(col,1.0),0.,1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *startunnel_source = "/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
"* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void lfnoise(in vec2 t, out float num);\n"
"void mfnoise(in vec2 x, in float fmin, in float fmax, in float alpha, out float num);\n"
"void dbox(in vec2 p, in vec2 b, out float dst);\n"
"void dpolygon(in vec2 x, in float N, in float R, out float dst);\n"
"void dstar(in vec2 x, in float N, in vec2 R, out float dst);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"void stroke(in float d0, in float s, out float d);\n"
"\n"
"// Mix appropriate marble colors.\n"
"void color(in float scale, out vec3 col)\n"
"{\n"
"    const int N = 13;\n"
"    const vec3 colors[N] = vec3[N](\n"
"        c.yyy,\n"
"        vec3(0.15,0.14,0.12),\n"
"        vec3(0.38,0.16,0.16),\n"
"        vec3(0.42,0.20,0.19),\n"
"        vec3(0.60,0.14,0.16),\n"
"        vec3(0.70,0.11,0.15),\n"
"        vec3(0.89,0.11,0.10),\n"
"        vec3(0.89,0.27,0.03),\n"
"        vec3(0.92,0.39,0.14),\n"
"        vec3(0.91,0.47,0.15),\n"
"        vec3(0.92,0.57,0.14),\n"
"        vec3(0.90,0.63,0.12),\n"
"        vec3(0.92,0.72,0.14)\n"
"    );\n"
"	float index = floor(scale*float(N)), \n"
"        remainder = scale*float(N)-index;\n"
"    col = mix(colors[int(index)],colors[int(index)+1], remainder);\n"
"}\n"
"\n"
"void colorize(in vec2 uv, out vec3 col, float i)\n"
"{\n"
"    vec2 n, n2;\n"
"    lfnoise(i*c.xx-2.*iTime, n.x);\n"
"    lfnoise(i*c.xx-2.*iTime-1337., n.y);\n"
"    \n"
"    uv += i*.5*n;\n"
"    \n"
"    float ca = .5*i-3.5*iTime,\n"
"        cc = cos(ca),\n"
"        sc = sin(ca);\n"
"    mat2 RR = mat2(cc,sc,-sc,cc);\n"
"    uv = RR*uv;\n"
"    \n"
"    float dd;\n"
"    rand(floor(.33*iTime)*c.xx,dd);\n"
"    \n"
"    vec3 c1;\n"
"    color(clamp(i+2.*n.y,0.,1.), c1);\n"
"    \n"
"    float d, da, db;\n"
"    dpolygon(uv, max(ceil(8.*dd),3.), .4+0.*.4*iScale, d);\n"
"    dstar(uv,max(ceil(8.*dd),3.),vec2(.05,.5)+0.*vec2(.1,.4)*iScale,db);\n"
"    d = mix(d,db,.5+.5*n.y);\n"
"    stroke(d, .01, da);\n"
"    da -= .01*n.y;\n"
"    d -= .01*n.x;\n"
"\n"
"        float mat = n.y;\n"
"//     rand(i*c.xx, mat);\n"
"        float phi = atan(uv.y, uv.x),\n"
"                dhex,\n"
"                na,\n"
"                nal;\n"
"            vec2 ind;\n"
"            rand(floor(.33*iTime)*c.xx, na);\n"
"            rand(floor(.33*iTime)*c.xx+1., nal);\n"
"            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));\n"
"            \n"
"            mat3 RRR;\n"
"            rot3(na*1.e3*vec3(1.1,1.5,1.9),RRR);\n"
"\n"
"            c1 = mix((.5+.5*n.y)*c1,(1.+.8*mat)*abs(RRR*c1),.5+.5*n.x);\n"
"            c1 = mix(c1,.2*c1,.5+.5*n.y);\n"
"\n"
"    \n"
"    col = mix(c1,c.yyy,smoothstep(1.5/iResolution.y, -1.5/iResolution.y,-da));\n"
"    col = mix(col, mix(col,c1,.03+.02*clamp(iScale,0.,1.)), smoothstep(1.5/iResolution.y, -1.5/iResolution.y,-d));\n"
"//     col = mix(col, mix(col,1.4*c1,.03+.02*clamp(iScale,0.,1.)), smoothstep(1.5/iResolution.y, -1.5/iResolution.y,abs(d)-.05));\n"
"    col = mix(col,(.5+iScale)*col,smoothstep(1.5/iResolution.y, -1.5/iResolution.y,-abs(da)+.005));\n"
"    \n"
"    col = clamp(col*1.3, 0.,1.);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    /// Set up coordinates\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(col, 0.);\n"
"        return;\n"
"    }\n"
"    \n"
"    vec3 x,\n"
"        o = c.yyx,\n"
"        t = vec3(uv,0.),\n"
"        dir = normalize(t-o),\n"
"        c1;\n"
"    float d = .5;\n"
"    vec2 n;\n"
"    int N = 100,\n"
"        i;\n"
"    \n"
"    for(i=0; i<N; ++i)\n"
"    {\n"
"        d = -(o.z-.5+.1*float(i))/dir.z;\n"
"        x = o + d * dir;\n"
"           \n"
"        colorize(x.xy,c1, float(i)*mix(.01,.05,iScale));\n"
"        col += c1;\n"
"    }\n"
"    c1 /= float(N);\n"
"    \n"
"    fragColor = vec4(col,1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *team210_logo_source = "/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
"* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void lfnoise(in vec2 t, out float num);\n"
"void mfnoise(in vec2 x, in float fmin, in float fmax, in float alpha, out float num);\n"
"\n"
"void camerasetup(in vec2 uv, out vec3 ro, out vec3 dir)\n"
"{\n"
"    vec3 right = c.xyy, up = c.yxy, target = c.yyy;\n"
"    ro = c.yyx+.3*vec3(cos(iTime), sin(iTime), 0.);\n"
"    dir = normalize(target + uv.x * right + uv.y * up - ro);\n"
"}\n"
"\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dcirclesegment(in vec2 x, in float R, in float p0, in float p1, out float d);\n"
"void dcircle(in vec2 x, in float R, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void dlogo210(in vec2 x, in float R, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void dhexagonpattern(in vec2 p, out float d, out vec2 ind);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"\n"
"// graph traversal for 210 logo effect\n"
"void textpre(in vec3 x, out vec2 sdf)\n"
"{\n"
"    float blend = 1.;\n"
"    dbox3(x, vec3(1., .5, .01+blend), sdf.x);\n"
"}\n"
"\n"
"// Perform raymarching for bounding object\n"
"void marchbounds(in vec3 ro, in vec3 dir, in int N, in float eps, out vec3 x, out vec2 s, out float d, out bool flag)\n"
"{\n"
"    flag = false;\n"
"    for(int ia=0; ia<max(N,0); ++ia)\n"
"	{\n"
"        x = ro + d*dir;\n"
"        textpre(x,s);\n"
"        if(s.x < eps)\n"
"        {\n"
"            flag = true;\n"
"            break;\n"
"        }\n"
"        d += s.x;\n"
"	}\n"
"}\n"
"\n"
"// 3D Effect on text in intro (210 logo)\n"
"void texteffect(in vec3 x, out vec2 sdf)\n"
"{\n"
"    // Start with z=0 plane\n"
"    sdf = vec2(x.z, 7.0);\n"
"    vec2 ind;\n"
"    float hex;\n"
"    dhexagonpattern(48.0*x.xy, hex, ind);\n"
"    \n"
"    // compute hexagon indices in cartesian coordinates\n"
"    vec2 cind = ind/48.0;\n"
"    \n"
"    // build up team210 logo (t < 12.)\n"
"    float inner_logo, logo_border; \n"
"    dlogo210(3.5*cind, 1., inner_logo);\n"
"    stroke(inner_logo, 0.35, inner_logo);\n"
"\n"
"    vec2 n;\n"
"    lfnoise(ind-2.*iTime, n.x);\n"
"    lfnoise(ind-2.*iTime-1337., n.y);\n"
"    \n"
"    // blend back to structure (t < 16., t > 12.)\n"
"    float blend = clamp(.5+.5*n.x, 0., 1.);\n"
"    \n"
"    if(inner_logo < 0.0 && blend >= 1.0e-3)\n"
"    {\n"
"        float noise;\n"
"        lfnoise(24.0*cind.xy-iTime, noise);\n"
"        zextrude(x.z,\n"
"                 1.5*x.z-inner_logo, \n"
"                 .5*(0.5+0.5*noise)*blend*(.1+.9*clamp(iScale,0.,1.)),\n"
"                 sdf.x);\n"
"        stroke(sdf.x, 0.05*blend, sdf.x);\n"
"        sdf.y = 7.0;\n"
"    }\n"
"    stroke(sdf.x,0.1,sdf.x);\n"
"    \n"
"//     // Add guard objects for debugging\n"
"//     float dr = .02;\n"
"//     vec3 y = mod(x,dr)-.5*dr;\n"
"//     float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));\n"
"//     guard = abs(guard)+dr*.1;\n"
"//     sdf.x = min(sdf.x, guard);\n"
"}\n"
"\n"
"// Perform raymarching for actual object\n"
"void marchscene(in vec3 ro, in vec3 dir, in int N, in float eps, out vec3 x, out vec2 s, out float d, out bool flag)\n"
"{\n"
"    flag = false;\n"
"    for(int ia=0; ia<max(N,0); ++ia)\n"
"	{\n"
"        x = ro + d*dir;\n"
"        texteffect(x,s);\n"
"        if(s.x < eps)\n"
"        {\n"
"            flag = true;\n"
"            break;\n"
"        }\n"
"        d += min(s.x,.005);\n"
"	}\n"
"}\n"
"\n"
"void calcnormal(in vec3 x, in float eps, out vec3 n)\n"
"{\n"
"    vec2 s, sp;\n"
"    texteffect(x, s);\n"
"    texteffect(x+eps*c.xyy, sp);\n"
"    n.x = sp.x-s.x;\n"
"    texteffect(x+eps*c.yxy, sp);\n"
"    n.y = sp.x-s.x;\n"
"    texteffect(x+eps*c.yyx, sp);\n"
"    n.z = sp.x-s.x;\n"
"    n = normalize(n);\n"
"}\n"
"\n"
"// Initial intro\n"
"vec2 ind;\n"
"void background2(in vec2 uv, out vec3 col)\n"
"{\n"
"    col = c.yyy;\n"
"    \n"
"    // hexagonal grid\n"
"    float d, d0;\n"
"//     vec2 ind;\n"
"    dhexagonpattern(48.0*uv, d0, ind);\n"
"    d = -d0;\n"
"    stroke(d, 0.1, d);\n"
"    vec2 cind = ind/48.0;\n"
"    \n"
"    // build up team210 logo (t < 12.)\n"
"    float inner_logo, logo_border; \n"
"    dlogo210(3.5*cind, 1., inner_logo);\n"
"    stroke(inner_logo, 0.35, inner_logo);\n"
"    stroke(inner_logo, 0.08, logo_border);\n"
"    \n"
"    vec2 n;\n"
"    lfnoise(ind-2.*iTime, n.x);\n"
"    lfnoise(ind-2.*iTime-1337., n.y);\n"
"    \n"
"    // blend back to structure (t < 16., t > 12.)\n"
"    float blend = 0.;\n"
"    \n"
"    inner_logo = mix(inner_logo, d0, blend);\n"
"    logo_border = mix(logo_border, d0, blend);\n"
"\n"
"    // make background change the color with time\n"
"    vec2 dt;\n"
"    lfnoise(15.0*cind+2.0, dt.x);\n"
"    lfnoise(15.0*cind+3.0, dt.y);\n"
"    dt *= 2.;\n"
"    float dm, dm2;\n"
"    lfnoise(50.0*cind, dm);\n"
"    dm = 0.5+0.5*dm;\n"
"    lfnoise(6.5*cind-dt-2.0*iTime*c.xx, dm2);\n"
"    dm2 = 0.5+0.5*dm2;\n"
"    \n"
"    // Colors\n"
"    vec3 orange = vec3(0.26,0.35,0.40);\n"
"    orange = mix(vec3(.60,0.24,0.60), orange, dm2);\n"
"    vec3 gray = .25*length(orange)/sqrt(3.)*c.xxx;\n"
"    \n"
"    float mat = n.x;\n"
"//     rand(i*c.xx, mat);\n"
"    float phi = atan(uv.y, uv.x),\n"
"        dhex,\n"
"        na,\n"
"        nal;\n"
"    rand(floor(.33*iTime)*c.xx, na);\n"
"    rand(floor(.33*iTime)*c.xx+1., nal);\n"
"    na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));\n"
"    \n"
"    mat3 RRR;\n"
"    rot3(na*1.e3*vec3(1.1,1.5,1.9)+iNote*210.,RRR);\n"
"\n"
"    orange = mix((.5+.5*n.y)*orange,(1.+.8*mat)*abs(RRR*orange),.5+.5*n.x);\n"
"    orange = mix(orange,.2*RRR*orange,.5+.5*n.y);\n"
"  \n"
"    col = orange;\n"
"    col = mix(col, gray, step(0.,inner_logo));\n"
"    \n"
"    col = mix(col, mix(gray,orange,step(abs(d0)-.4,0.)), step(-abs(d0)+.2,0.));\n"
"    col = mix(col, 4.*orange, step(logo_border,0.));\n"
"    \n"
"    col = 2.*col;\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), s = c.xy;\n"
"\n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(c.yyy, 0.);\n"
"        return;\n"
"    }\n"
"    uv *= 1.7;\n"
"    \n"
"	vec3 ro, x, dir;\n"
"    \n"
"    float d = 0.;\n"
"    bool hit = false;\n"
"    \n"
"    vec3 col = c.yyy;\n"
"                \n"
"	camerasetup(uv, ro, dir);\n"
"    d = (.5-ro.z)/dir.z;\n"
"    marchbounds(ro, dir, 150, 2.0e-4, x, s, d, hit);\n"
"\n"
"    if(hit) hit = false;\n"
"    else d = -ro.z/dir.z;\n"
"    marchscene(ro, dir, 800, 2.0e-4, x, s, d, hit);\n"
"    \n"
"    if(hit)\n"
"    {\n"
"        vec3 n;\n"
"        calcnormal(x, 2.0e-4, n);\n"
"\n"
"        float rs = 1.9;\n"
"        vec3 l = x+1.*c.yyx,\n"
"        	re = normalize(reflect(-l,n));\n"
"        float rev = abs(dot(re,dir)), ln = abs(dot(l,n));\n"
"		background2(x.xy, col);\n"
"        col = mix(col, .5*col, clamp(x.z/.2,0.,1.));\n"
"    }\n"
"    else\n"
"    {\n"
"        float ra;\n"
"        rand(ind, ra);\n"
"        background2((ro-(ro.z/dir.z-.5*ra)*dir).xy, col);\n"
"    }\n"
"    col = clamp(col, 0., 1.);\n"
"    fragColor = vec4(col,1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *broccoli_source = "/* Corfield Imitation 1\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" * \n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" * \n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" * \n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
" */\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void lfnoise(in vec2 t, out float n);\n"
"\n"
"void rand3(in vec3 x, out float num);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"\n"
"// Random Quadtree\n"
"void dcubetree(in vec3 x, in float threshold, in float depth, out float d, out float faco)\n"
"{\n"
"    d = 1.;\n"
"    vec3 y = x, \n"
"        yi;\n"
"    float size = .5,\n"
"	    fac = 1.;\n"
"    faco = 1.;\n"
"    for(float i=0.; i<depth; i+=1.)\n"
"    {\n"
"        vec3 y0 = y;\n"
"        y = mod(y, size)-.5*size;\n"
"        yi = y0-y;\n"
"		float r;\n"
"        rand3(yi+fac,r);\n"
"        fac *= r*step(r,threshold);\n"
"        if(fac != 0.)\n"
"        {\n"
"            float dd;\n"
"            dbox3(y,(.35)*size*c.xxx,dd);\n"
"//             dd = mix(dd,length(y)-(.3)*size,step(r,threshold));\n"
"            dd = abs(dd)-.01*size;\n"
"            smoothmin(d,dd,.01,d);\n"
"        } else break;\n"
"        \n"
"        size *= .5;\n"
"    }\n"
"    faco += fac*fac;\n"
"}\n"
"\n"
"void rot3(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\n"
"// Scene\n"
"float mat;\n"
"void scene(in vec3 x, out vec2 d)\n"
"{\n"
"    d = c.xx;\n"
"    \n"
"    x.z -= mix(.01,.1,iDial0)*iTime;\n"
"    \n"
"    dcubetree((2.)*x-iTime*c.yyx-.1*iTime, .5,  6.-6.*iScale, d.x, mat);\n"
"    float d2, m2 = 0.;\n"
"    //dcubetree((3.+.5*iScale)*x+iTime*c.yyx-.1*iTime, .5, 6., d2, m2);\n"
"//     mat = 5555.*mat;\n"
"	//lfnoise((5.*x.z-5.*iTime)*c.xx, mat);\n"
"	//mat = .5+.5*mat;\n"
"    //smoothmin(d.x,d2,.2+.2*iScale,d.x);\n"
"    \n"
"    d2 = length(x.xy)-.1;\n"
"    d.x = max(d.x, -d2);\n"
"    d = min(d, -length(x.xy)+.3);\n"
"    \n"
"//     d -= .01;\n"
"    \n"
"    stroke(d.x, .01, d.x);\n"
"}\n"
"\n"
"// Normal\n"
"const float dx = 5.e-4;\n"
"void normal(in vec3 x, out vec3 n);\n"
"\n"
"// Texture\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{    \n"
"    float phi = .1*iTime;\n"
"    \n"
"    vec3 white = vec3(0.89,0.44,0.23),\n"
"        gray =vec3(0.25,0.23,0.21);\n"
"    float size = .1;\n"
"    \n"
"    \n"
"    vec2 y = mod(x,size)-.5*size;\n"
"    y = abs(y)-.001;\n"
"    \n"
"    float d = min(y.x,y.y);\n"
"    col = mix(white, gray, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d*mat));\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    // Set up coordinates\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(col, 0.);\n"
"        return;\n"
"    }\n"
"        \n"
"     // Camera setup\n"
"    float pp = .3*iTime;\n"
"    vec3 o = c.yyx,\n"
"        t = c.yyy,\n"
"        dir = normalize(t-o),\n"
"        r = normalize(c.xyy),\n"
"        u = normalize(cross(r,dir)),\n"
"        n,\n"
"        x,\n"
"        l;\n"
"    t += uv.x*r + uv.y*u;\n"
"    dir = normalize(t-o);\n"
"    vec2 s;\n"
"    float d = (.1)/length(dir.xy);// -(o.z-.12)/dir.z;\n"
"    int N = 650,\n"
"        i;\n"
"    \n"
"    // Graph\n"
"    x = o + d * dir;\n"
"    \n"
"    // Actual Scene\n"
"    {\n"
"\n"
"        // Raymarching\n"
"        for(i=0; i<N; ++i)\n"
"        {\n"
"            x = o + d * dir;\n"
"            scene(x,s);\n"
"            if(s.x < 1.e-4) break;\n"
"            d += min(s.x,.0005);\n"
"        }\n"
"\n"
"        // Illumination\n"
"        l = normalize(x+c.yxx);\n"
"        if(i<N)\n"
"        {\n"
"            normal(x,n);\n"
"            mat += 14.*sign(n.z);\n"
"            col = mix((.5+.5*mat)*c.xxx,(1.+.8*mat)*vec3(0.89,0.44,0.23),.5+.5*sin(x.z));\n"
"            col = mix(col,vec3(0.25,0.23,0.21),.5+.5*cos(4.*x.z+mat));\n"
"            float phi = atan(x.y, x.x),\n"
"                dhex,\n"
"                na,\n"
"                nal;\n"
"            vec2 ind;\n"
"            rand(floor(.33*iTime)*c.xx, na);\n"
"            rand(floor(.33*iTime)*c.xx+1., nal);\n"
"            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));\n"
"            \n"
"            mat3 RR;\n"
"            float ras;\n"
"            rand(mat*c.xx,ras);\n"
"            rot3(na*1.e3*vec3(1.1,1.5,1.9)+1.*ras+.5*cos(x.z)+iNote*210.,RR);\n"
"\n"
"            col = mix((.5+.5*mat)*c.xxx,(1.+.8*mat)*abs(RR*vec3(0.89,0.44,0.23)),.5+.5*sin(x.z));\n"
"            //rot3(c.xxx+x.z+1200.*na*mat+1.e3*na,RR);\n"
"            col = mix(col,abs(RR*vec3(0.25,0.23,0.21)),.5+.5*cos(.5*(x.z)));\n"
"            \n"
"            col = mix(col, .5*abs(RR*RR*vec3(0.25,0.23,0.21)), clamp(length(x.xy)/.2, 0.,1.));\n"
"            col = mix(col, abs(RR*col), step(0.,length(sign(n))));\n"
"            \n"
"            col = mix(col, 3.*col, step(.9,length(x.xy))*step(1.1,length(x.xy)));\n"
"            col = mix(col, .3*length(col)/sqrt(3.)*c.xxx, iPressure);\n"
"        }\n"
"    }\n"
"    vec3 c1 = col;\n"
"    // Colorize\n"
"    col = .8*col\n"
"        + .9*col*abs(dot(l,n))\n"
"        +5.4*col*abs(pow(dot(reflect(-l,n),dir),3.));\n"
"        \n"
"    float dd;\n"
"    rand(1200.*uv, dd);\n"
"    col += dd*.1*c.xxx;\n"
"    \n"
"    col = mix(col, c1, clamp(float(i)/float(N),0.,1.));\n"
"    \n"
"    fragColor = clamp(vec4(col,1.0),0.,1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *boxplane_source = "#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void rand(in vec2 x, out float n)\n"
"{\n"
"    x += 400.;\n"
"    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\n"
"void lfnoise(in vec2 t, out float n)\n"
"{\n"
"    vec2 i = floor(t);\n"
"    t = fract(t);\n"
"    t = smoothstep(c.yy, c.xx, t);\n"
"    vec2 v1, v2;\n"
"    rand(i, v1.x);\n"
"    rand(i+c.xy, v1.y);\n"
"    rand(i+c.yx, v2.x);\n"
"    rand(i+c.xx, v2.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    n = mix(v1.x, v1.y, t.x);\n"
"}\n"
"\n"
"void rot3(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\n"
"// Box sdf\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"    vec2 da = abs(x)-b;\n"
"    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"/*\n"
"//distance to spline with parameter t\n"
"float dist2(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t)\n"
"{\n"
"    t = clamp(t, 0., 1.);\n"
"    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);\n"
"}\n"
"\n"
"//minimum dist3ance to spline\n"
"void dspline2(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float ds)\n"
"{\n"
"    //coefficients for 0 = t^3 + a * t^2 + b * t + c\n"
"    vec2 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;\n"
"    vec3 ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);\n"
"\n"
"	//discriminant and helpers\n"
"    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;\n"
"    \n"
"    //triple real root\n"
"    if(dis > 0.) \n"
"    {\n"
"        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);\n"
"        ds = dist2(p0,p1,p2,x,ui.x+ui.y-tau);\n"
"        return;\n"
"    }\n"
"    \n"
"    //three dist3inct real roots\n"
"    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;\n"
"    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;\n"
"    ds = min(\n"
"        dist2(p0,p1,p2,x, t.x),\n"
"        min(\n"
"            dist2(p0,p1,p2,x,t.y),\n"
"            dist2(p0,p1,p2,x,t.z)\n"
"        )\n"
"    );\n"
"}\n"
"\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}*/\n"
"\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)\n"
"{\n"
"    vec3 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\n"
"// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\n"
"// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\n"
"// Add sdfs\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf)\n"
"{\n"
"    sdf = mix(sda, sdb, step(sdb.x, sda.x));\n"
"}\n"
"\n"
"vec2 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.y += mix(.3,5.,iDial0)*iTime;\n"
"    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));\n"
"    x.xy = R*x.xy;\n"
"    \n"
"    float d,\n"
"        size = mix(.5,.01,iDial6);\n"
"    vec2 x2 = mod(x.xy,size)-.5*size;\n"
"	\n"
"    ind = (x.xy - x2)/size;\n"
"    dbox(x2, .5*size*c.xx, d);\n"
"    zextrude(x.z, -d-.005, .05, d);\n"
"    d = max(x.z,d);\n"
"    d = abs(d);\n"
"    sdf = vec2(d,2.);\n"
"    \n"
"    float r;\n"
"    rand(ind-1.e2*floor(iTime), r);\n"
"    //lfnoise(12.*ind-1.*iTime, r);\n"
"    //r = .5+.5*r;\n"
"    if(r > .7)\n"
"    {\n"
"        dbox(x2, .5*size*c.xx, d);\n"
"        zextrude(x.z, -d-.02, (.1+.15*iScale)*(r-.7)/.3, d);\n"
"        stroke(d, .001, d);\n"
"        add(sdf, vec2(d,1.), sdf);\n"
"    }\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n)\n"
"{\n"
"    const float dx = 1.e-4;\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    x.y += mix(.3,5.,iDial0)*iTime;\n"
"    mat2 R = mat2(cos(pi/4.), sin(pi/4.), -sin(pi/4.), cos(pi/4.));\n"
"    x.xy = R*x.xy;\n"
"    \n"
"    float d,\n"
"        size = mix(.5,.01,iDial6),\n"
"        r;\n"
"    vec2 x2 = mod(x.xy,size)-.5*size;\n"
"    \n"
"    rand(ind-1.e2*floor(iTime), r);\n"
"    //lfnoise(12.*ind-1.*iTime, r);\n"
"    //r = .5+.5*r;\n"
"    col = mix(.14*c.xxx, .33*c.xxx, r);\n"
"    dbox(x2, .35*size*c.xx, d);\n"
"    if(r > .9)\n"
"    {\n"
"        col = mix(col, mix(c.xxy, c.xxx, .8), sm(d));\n"
"        stroke(d, .0025, d);\n"
"        col = mix(col, mix(c.xyy,c.xxx,.8), sm(d));\n"
"        stroke(d-.004, .002, d);\n"
"        col = mix(col, c.xyy, sm(d));\n"
"    }\n"
"	else if(r > .8)\n"
"    {\n"
"        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));\n"
"        stroke(d, .0025, d);\n"
"        col = mix(col, mix(.7*c.xxy,c.xxx,.8), sm(d));\n"
"        stroke(d-.004, .002, d);\n"
"        col = mix(col, .7*c.xxy, sm(d));\n"
"    }\n"
"    else if(r > .7)\n"
"    {\n"
"        col = mix(col, mix(c.xyy, c.xxx, .8), sm(d));\n"
"        stroke(d, .0025, d);\n"
"        col = mix(col, mix(mix(c.xxy, c.xyy, .5),c.xxx,.8), sm(d));\n"
"        stroke(d-.004, .002, d);\n"
"        col = mix(col, mix(c.xxy, c.xyy, .5), sm(d));\n"
"    }\n"
"    \n"
"    // Truchet\n"
"    /*\n"
"    float da = floor(4.*r)*pi/2.;\n"
"    R = mat2(cos(da), sin(da),-sin(da),cos(da));\n"
"    x2 = R * x2;\n"
"    if(r > .3)\n"
"    {\n"
"    	dspline2(x2,-.5*size*c.xy, c.yy, -.5*size*c.yx, d);\n"
"        dspline2(x2,.5*size*c.xy, c.yy, .5*size*c.yx, da);\n"
"        \n"
"    }\n"
"    else\n"
"    {\n"
"        dlinesegment(x2,-.5*size*c.xy, .5*size*c.xy, d);\n"
"        dlinesegment(x2,-.5*size*c.yx, .5*size*c.yx, da);\n"
"    }\n"
"    d = min(d, da);\n"
"    stroke(d, .001, d);\n"
"    col = mix(col,  .9*c.xyy, sm(d));\n"
"    stroke(d-.004, .002, d);\n"
"    col = mix(col, .0*c.xyy, sm(d));\n"
"*/\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx), \n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 100,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = -(o.z-.15)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(x.z<-.05)\n"
"        {\n"
"            col = .2*c.xxx;\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += min(s.x,5.e-3);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            vec3 c1;\n"
"            \n"
"            float r;\n"
"		    rand(ind-1.e2*floor(iTime), r);\n"
"            if(r > .9)\n"
"                col = c.xyy;\n"
"            else if(r > .8)\n"
"                col = .7*c.xxy;\n"
"            else if(r > .7)\n"
"                col = mix(c.xxy, c.xyy, .5);\n"
"            float sc = clamp((r-.7)/.3,0.,1.);\n"
"            col = mix(mix(col, c.xxx, .1*sc), .4*c.xyy, sc);\n"
"            col = .3*col\n"
"                + .9*col * abs(dot(l,n))\n"
"                + 1.3*col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"            col = mix(col, c.xxx, .4);\n"
"            col *= col;\n"
"            \n"
"            d = -(o.z)/dir.z;\n"
"            x = o + d * dir;\n"
"            scene(x,s);\n"
"            l = normalize(x+c.xzx);\n"
"            colorize(x.xy, c1);\n"
"            n = c.yyx;\n"
"            \n"
"            c1 = .1*c1\n"
"                + .8*c1 * abs(dot(l,n))\n"
"                + c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"            col = mix(col, c1, .3);\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            float r;\n"
"            \n"
"            colorize(x.xy, col);\n"
"            col = .1*col\n"
"                + .8*col * abs(dot(l,n))\n"
"                + col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"        }\n"
"    }\n"
"    col += col;\n"
"    col *= col;\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *doublependulum_source = "/* Corfield Imitation 1\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" * \n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" * \n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" * \n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
" */\n"
"  \n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"//lengths\n"
"float L1 = 2.;\n"
"float L2 = 1.;\n"
"\n"
"//masses\n"
"float M1 = 3.5;\n"
"float M2 = 0.1;\n"
"\n"
"//constants\n"
"float G = 9.81;\n"
"\n"
"//runge kutta params\n"
"float h = 1.e-1;\n"
"float tmax = 3.;\n"
"\n"
"/**eval system of differential equations\n"
"params:\n"
"tp[0]: theta1\n"
"tp[1]: theta2\n"
"tp[2]: ptheta1\n"
"tp[3]: ptheta2\n"
"*/\n"
"vec4 f(vec4 tp)\n"
"{\n"
"    float p0w = (L1*L2*(M1+M2*pow(sin(tp[0]-tp[1]),2.)));\n"
"    float C1 = tp[2]*tp[3]*sin(tp[0]-tp[1])/p0w;\n"
"    float C2 = (L2*L2*M2*tp[2]*tp[2]+L1*L1*(M1+M2)*tp[3]*tp[3]-L1*L2*M2*tp[2]*tp[3]*cos(tp[0]-tp[1]))*sin(2.*(tp[0]-tp[1]))/(2.*p0w*p0w);\n"
"    \n"
"    vec4 ret;\n"
"    \n"
"    ret[0] = (L2*tp[2]-L1*tp[3]*cos(tp[0]-tp[1]))/(L1*p0w);\n"
"    ret[1] = (L1*(M1+M2)*tp[3]-L2*M2*tp[2]*cos(tp[0]-tp[1]))/(L2*M2*p0w);\n"
"    ret[2] = -(M1+M2)*G*L1*sin(tp[0])-C1+C2;\n"
"    ret[3] = -M2*G*L2*sin(tp[1])+C1-C2;\n"
"    \n"
"    return ret;    \n"
"}\n"
"\n"
"vec4 step_rk4(vec4 tp)\n"
"{\n"
"    vec4 k1 = f(tp);\n"
"    vec4 k2 = f(tp + h/2.*k1);\n"
"    vec4 k3 = f(tp + h/2.*k2);\n"
"    vec4 k4 = f(tp + h*k3);\n"
"    return tp + h/6.*(k1+2.*k2+2.*k3+k4);\n"
"}\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"\n"
"// Scene\n"
"float mat;\n"
"void scene(in vec3 x, out vec2 d)\n"
"{\n"
"    d = c.xx;\n"
"    \n"
"    x.y -= .1*iTime;\n"
"    x.x += mix(0.,10., iDial6);\n"
"    \n"
"    vec4 state = vec4(x.xy*2.*pi-vec2(pi,pi), 0, 0);\n"
"    float time = 0.;\n"
"    while(time < tmax) \n"
"    {\n"
"        state = step_rk4(state);\n"
"        time += h;\n"
"    }\n"
"    \n"
"    d.x = x.z -.02 + .01*abs(log(abs(state.r/state.b)));\n"
"    //stroke(d.x, .005, d.x);\n"
"}\n"
"\n"
"// Normal\n"
"const float dx = 5.e-4;\n"
"void normal(in vec3 x, out vec3 n);\n"
"\n"
"// Texture\n"
"void color(in float scale, out vec3 col)\n"
"{\n"
"    scale = fract(scale);\n"
"    const int N = 5;\n"
"    const vec3 colors[N] = vec3[N](\n"
"        vec3(1.00,0.67,0.36),\n"
"        vec3(0.86,0.45,0.50),\n"
"        vec3(0.67,0.43,0.51),\n"
"        vec3(0.41,0.36,0.48),\n"
"        vec3(0.27,0.36,0.48)\n"
"    );\n"
"	float index = floor(scale*float(N)), \n"
"        remainder = scale*float(N)-index;\n"
"    col = mix(colors[int(index)],colors[int(index)+1], remainder);\n"
"    col = fract(col);\n"
"}\n"
"\n"
"// Texture\n"
"void color2(in float scale, out vec3 col)\n"
"{\n"
"    scale = fract(scale);\n"
"    const int N = 5;\n"
"    const vec3 colors[N] = vec3[N](\n"
"        vec3(0.95,0.74,0.56),\n"
"        vec3(0.95,0.64,0.50),\n"
"        vec3(0.75,0.43,0.42),\n"
"        vec3(0.44,0.29,0.36),\n"
"        vec3(0.24,0.15,0.23)\n"
"    );\n"
"	float index = floor(scale*float(N)), \n"
"        remainder = scale*float(N)-index;\n"
"    col = mix(colors[int(index)],colors[int(index)+1], remainder);\n"
"    col = fract(col);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    // Set up coordinates\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"\n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(col, 0.);\n"
"        return;\n"
"    }\n"
"    vec3 na;\n"
"    lfnoise(iTime*c.xx, na.x);\n"
"    lfnoise(2.*iTime*c.xx+2337., na.y);\n"
"    lfnoise(3.*iTime*c.xx+3337., na.z);\n"
"    \n"
"    M1 = 2.+mix(0.,1.,iScale)*na.x;\n"
"    M2 = 2.+mix(0.,1.,iScale)*na.y;\n"
"    \n"
"    L1 = 2.+.1*mix(0.,1.,iScale)*na.x;\n"
"    L2 = 1.+.1*mix(0.,1.,iScale)*na.z;\n"
"//     h = mix(1.e-1,1., iDial7);\n"
"    tmax = mix(1.,4., iDial7);\n"
"    uv /= mix(.5,8.,iDial0);\n"
"    \n"
"     // Camera setup\n"
"    float pp = .3*iTime;\n"
"    vec3 o = c.yyx,\n"
"        t = c.yyy,\n"
"        dir = normalize(t-o),\n"
"        r = normalize(c.xyy),\n"
"        u = normalize(cross(r,dir)),\n"
"        n,\n"
"        x,\n"
"        l;\n"
"    t += uv.x*r + uv.y*u;\n"
"    dir = normalize(t-o);\n"
"    vec2 s;\n"
"    float d = -(o.z)/dir.z;\n"
"    int N = 450,\n"
"        i;\n"
"    \n"
"    // Graph\n"
"    x = o + d * dir;\n"
"    \n"
"    // Actual Scene\n"
"    {\n"
"\n"
"        // Raymarching\n"
"//         for(i=0; i<N; ++i)\n"
"//         {\n"
"//             x = o + d * dir;\n"
"//             scene(x,s);\n"
"//             if(s.x < 1.e-4) break;\n"
"//             d += s.x;//,.005);\n"
"//         }\n"
"\n"
"        // Illumination\n"
"        l = normalize(x+c.yxx);\n"
"        if(i<N)\n"
"        {\n"
"            vec4 state = vec4((x.xy-.1*iTime*c.yx+mix(0.,10., iDial6)*c.xy)*2.*pi-vec2(pi,pi), 0, 0);\n"
"            float time = 0.;\n"
"            while(time < tmax) \n"
"            {\n"
"                state = step_rk4(state);\n"
"                time += h;\n"
"            }\n"
"            float da = -.02 - mix(.01,.03,iScale)*log(abs(state.r/state.b));\n"
"            d += da;\n"
"            x = o + d * dir;\n"
"            normal(x,n);\n"
"            color(abs(x.z-.1)*1., col);\n"
"            vec3 c1 = c.yyy;\n"
"            color2(abs(x.z-.1)*1., c1);\n"
"            \n"
"            col = mix(col, c1, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, abs(da)-.1));\n"
"            \n"
"                        float na, nal;\n"
"			rand(floor(.33*iTime)*c.xx, na);\n"
"            rand(floor(.33*iTime)*c.xx+1., nal);\n"
"            na = mix(na,nal,clamp(((.33*iTime-floor(.33*iTime))-.9)/.1,0.,1.));\n"
"             mat3 RR;\n"
"            rot3(na*1.e3*vec3(1.1,1.5,1.9)+13.*length(col)+iNote*210.,RR);\n"
"            col = abs(RR*col);\n"
"            \n"
"            col = mix(col, length(col)/sqrt(3)*c.xxx, .7);\n"
"        }\n"
"    }\n"
"    \n"
"    // Colorize\n"
"    col = .2*col\n"
"        + 1.3*col*abs(dot(l,n))\n"
"        +.4*col*abs(pow(dot(reflect(-l,n),dir),3.));\n"
"    \n"
"    float dd;\n"
"    rand(1200.*uv, dd);\n"
"    col += dd*.1*c.xxx;\n"
"    \n"
"    fragColor = clamp(vec4(col,1.0),0.,1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *midikeyboard_source = "#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"    vec2 da = abs(x)-b;\n"
"    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\n"
"void key(in float keycode, out float windex, out float bindex)\n"
"{\n"
"    float note = mod(keycode-48., 12.),\n"
"        oct = round((keycode-48.-note)/12.);\n"
"    if(note == 0.) //C\n"
"    {\n"
"        windex = 7.*oct;\n"
"        bindex = -499.;\n"
"    }\n"
"    else if(note == 1.) //C#\n"
"    {\n"
"        windex = -499.;\n"
"        bindex = 7.*oct;\n"
"    }\n"
"    else if(note == 2.) //D\n"
"    {\n"
"        windex = 7.*oct + 1.;\n"
"        bindex = -499.;\n"
"    }\n"
"	else if(note == 3.) //D#\n"
"    {\n"
"        windex = -499.;\n"
"        bindex = 7.*oct + 1.;\n"
"    }\n"
"    else if(note == 4.) //E\n"
"    {\n"
"        windex = 7.*oct + 2.;\n"
"        bindex = -499.;\n"
"    }\n"
"	else if(note == 5.) //F\n"
"    {\n"
"        windex = 7.*oct + 3.;\n"
"        bindex = -499.;\n"
"    }\n"
"   	else if(note == 6.) //F#\n"
"    {\n"
"        windex = -499.;\n"
"        bindex = 7.*oct + 3.;\n"
"    }\n"
"    else if(note == 7.) //G\n"
"    {\n"
"        windex = 7.*oct + 4.;\n"
"        bindex = -499.;\n"
"    }\n"
"   	else if(note == 8.) //G#\n"
"    {\n"
"        windex = -499.;\n"
"        bindex = 7.*oct + 4.;\n"
"    }\n"
"    else if(note == 9.) //A\n"
"    {\n"
"        windex = 7.*oct + 5.;\n"
"        bindex = -499.;\n"
"    }\n"
"   	else if(note == 10.) //A#\n"
"    {\n"
"        windex = -499.;\n"
"        bindex = 7.*oct + 5.;\n"
"    }\n"
"    else if(note == 11.) //H\n"
"    {\n"
"        windex = 7.*oct + 6.;\n"
"        bindex = -499.;\n"
"    }\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 x = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    x.x += 7.*a/16.;\n"
"    \n"
"//     iNote = 62.;\n"
"    \n"
"    // determine note index\n"
"    float windex, bindex;\n"
"    key(iNote, windex, bindex);\n"
"    \n"
"    float d;\n"
"    \n"
"    // white keys\n"
"    vec2 y = vec2(mod(x.x, a/16.)-a/32.,x.y);\n"
"    dbox(y-.2*c.yx, vec2(a/34.,.4), d);\n"
"    float inda = ceil((x.x-y.x)*16./a);\n"
"    if(inda < 24.)\n"
"	    col = mix(col, mix(c.xxx,c.yxy,float(inda==(windex))), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));\n"
"\n"
"    // black keys\n"
"    vec2 z = vec2(mod(x.x-a/32., a/16.)-a/32.,x.y);\n"
"    float ind = round((x.x-z.x)*16./a),\n"
"        cond = mod(ind, 7.);\n"
"	if(cond != 2. && cond != 6.)\n"
"    {\n"
"    	dbox(z-.4*c.yx-.0*c.xy, vec2(a/68.,.4), d);\n"
"    	col = mix(col, mix(c.yyy, c.yyx, float(ind == bindex)), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));\n"
"    }\n"
"    \n"
"    fragColor = vec4(col,1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *wursttunnel_source = "/* Corfield Imitation 1\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" * \n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" * \n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" * \n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
" */\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFFTWidth;\n"
"uniform float iScale;\n"
"uniform float iHighScale;\n"
"uniform float iNBeats;\n"
"uniform float iDial0;\n"
"uniform float iDial6;\n"
"uniform float iDial7;\n"
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
"uniform float iNote;\n"
"uniform float iPressure;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"\n"
"void dtruchet3(in vec3 x, in float size, out float d)\n"
"{\n"
"    vec3 y = mod(x, size)-.5*size,\n"
"        ind = (x-y)/size;\n"
"    \n"
"    vec3 r;\n"
"    lfnoise(ind.xy-.1*iTime, r.x);\n"
"    lfnoise(ind.yz-.1*iTime, r.y);\n"
"    lfnoise(ind.zx-.1*iTime, r.z);\n"
"    \n"
"    mat3 RR;\n"
"    rot3(floor(12.*r)*pi/2., RR);\n"
"    \n"
"    dspline3(RR*y, .5*size*c.zyy, c.yyy, .5*size*c.yzy, d);\n"
"    float da;\n"
"    dspline3(RR*y, .5*size*c.yyz, c.yyy, .5*size*c.xyy, da);\n"
"    d = min(d, da);\n"
"    dspline3(RR*y, .5*size*c.yxy, c.yyy, .5*size*c.yyx, da);\n"
"    d = min(d, da);\n"
"    \n"
"    vec3 na, nb;\n"
"    lfnoise(62.*x.x*c.xx-1.*iTime, na.x);\n"
"    lfnoise(72.*x.y*c.xx-1.3*iTime, na.y);\n"
"    lfnoise(-62.*x.x*c.xx+1.*iTime, nb.x);\n"
"    lfnoise(-72.*x.y*c.xx+1.3*iTime, nb.y);\n"
"    \n"
"    stroke(d, .125*size+mix(.005,.007,iScale)*na.x*na.y-mix(.005,01,iScale)*nb.x*nb.y, d);\n"
"}\n"
"\n"
"// Scene\n"
"float mat;\n"
"void scene(in vec3 x, out vec2 d)\n"
"{\n"
"	x.z -= mix(.1,1.,iDial0)*iTime;\n"
"    float n;\n"
"    lfnoise(3.*x.z*c.xx-.1*iTime, n);\n"
"    n = .5+.5*n;\n"
"    \n"
"    vec2 cs = vec2(cos(iTime+3.*n), sin(iTime+3.*n));\n"
"    mat2 r2 = mat2(cs.x,cs.y,-cs.y,cs.x);\n"
"    x.xy = r2 * x.xy;\n"
"    \n"
"    d = c.xx;\n"
"    dtruchet3(x, .1, d.x);\n"
"}\n"
"\n"
"// Normal\n"
"const float dx = 5.e-4;\n"
"void normal(in vec3 x, out vec3 n);\n"
"\n"
"// Texture\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{    \n"
"    float phi = .1*iTime;\n"
"    \n"
"    vec3 white = vec3(0.89,0.44,0.23),\n"
"        gray =vec3(0.25,0.23,0.21);\n"
"    float size = .1;\n"
"    \n"
"    \n"
"    vec2 y = mod(x,size)-.5*size;\n"
"    y = abs(y)-.001;\n"
"    \n"
"    float d = min(y.x,y.y);\n"
"    col = mix(white, gray, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d*mat));\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    // Set up coordinates\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(col, 0.);\n"
"        return;\n"
"    }\n"
"    \n"
"     // Camera setup\n"
"    float pp = .3*iTime;\n"
"    vec3 o = c.yyx,\n"
"        t = c.yyy,\n"
"        dir = normalize(t-o),\n"
"        r = normalize(c.xyy),\n"
"        u = normalize(cross(r,dir)),\n"
"        n,\n"
"        x,\n"
"        l;\n"
"    t += uv.x*r + uv.y*u;\n"
"    dir = normalize(t-o);\n"
"    vec2 s;\n"
"    float d = 0.;//-(o.z-.03)/dir.z;\n"
"    int N = 450,\n"
"        i;\n"
"    \n"
"    // Graph\n"
"    x = o + d * dir;\n"
"    mat3 RR;\n"
"            vec3 ra;\n"
"            rand(iNote*c.xx+3337., ra.x);\n"
"            rand(iNote*c.xx+1337., ra.y);\n"
"            rand(iNote*c.xx+2337., ra.z);\n"
"            rot3((iNote-48.)*810.*ra,RR);\n"
"    // Actual Scene\n"
"    {\n"
"\n"
"        // Raymarching\n"
"        for(i=0; i<N; ++i)\n"
"        {\n"
"            x = o + d * dir;\n"
"            scene(x,s);\n"
"            if(s.x < 1.e-4) break;\n"
"            d += min(s.x,.01);\n"
"        }\n"
"\n"
"        // Illumination\n"
"        l = normalize(x-.1*c.yyx);\n"
"        if(i<N)\n"
"        {\n"
"            normal(x,n);\n"
"            \n"
"            vec3 lo = vec3(0.69,0.23,0.22),\n"
"                hi = vec3(1.00,0.33,0.27);\n"
"            col  = mix(hi, lo, tanh(d/3.));\n"
"            \n"
"            col = abs(RR*col);\n"
"            col = mix(col, .3*c.xxx, .6);\n"
"//             col = mix(col, .3*length(col)/sqrt(3.)*c.xxx, iPressure);\n"
"//             col = mix(col, c.yyy, clamp(float(i)/float(N),0.,1.));\n"
"        }\n"
"    }\n"
"    \n"
"    // Colorize\n"
"    col = .2*col\n"
"        + 1.8*col*abs(dot(l,n))\n"
"        +1.6*abs(RR*vec3(1.00,0.33,0.27))*(pow(dot(reflect(-l,n),dir),9.));\n"
"    col = clamp(col, 0.,1.);\n"
"    if(length(col)<.01)col = .1*abs(RR*vec3(0.69,0.23,0.22));\n"
"        \n"
"    float dd;\n"
"    rand(1200.*uv, dd);\n"
"    col += dd*.1*c.xxx;\n"
"    \n"
"    fragColor = clamp(vec4(col,1.0),0.,1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
void Loadrand()
{
    int rand_size = strlen(rand_source);
    rand_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rand_handle, 1, (GLchar **)&rand_source, &rand_size);
    glCompileShader(rand_handle);
#ifdef DEBUG
    printf("---> rand Shader:\n");
    debug(rand_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadzextrude()
{
    int zextrude_size = strlen(zextrude_source);
    zextrude_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(zextrude_handle, 1, (GLchar **)&zextrude_source, &zextrude_size);
    glCompileShader(zextrude_handle);
#ifdef DEBUG
    printf("---> zextrude Shader:\n");
    debug(zextrude_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadstroke()
{
    int stroke_size = strlen(stroke_source);
    stroke_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(stroke_handle, 1, (GLchar **)&stroke_source, &stroke_size);
    glCompileShader(stroke_handle);
#ifdef DEBUG
    printf("---> stroke Shader:\n");
    debug(stroke_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadsmoothmin()
{
    int smoothmin_size = strlen(smoothmin_source);
    smoothmin_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(smoothmin_handle, 1, (GLchar **)&smoothmin_source, &smoothmin_size);
    glCompileShader(smoothmin_handle);
#ifdef DEBUG
    printf("---> smoothmin Shader:\n");
    debug(smoothmin_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddhexagonpattern()
{
    int dhexagonpattern_size = strlen(dhexagonpattern_source);
    dhexagonpattern_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dhexagonpattern_handle, 1, (GLchar **)&dhexagonpattern_source, &dhexagonpattern_size);
    glCompileShader(dhexagonpattern_handle);
#ifdef DEBUG
    printf("---> dhexagonpattern Shader:\n");
    debug(dhexagonpattern_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadnormal()
{
    int normal_size = strlen(normal_source);
    normal_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(normal_handle, 1, (GLchar **)&normal_source, &normal_size);
    glCompileShader(normal_handle);
#ifdef DEBUG
    printf("---> normal Shader:\n");
    debug(normal_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrot3()
{
    int rot3_size = strlen(rot3_source);
    rot3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rot3_handle, 1, (GLchar **)&rot3_source, &rot3_size);
    glCompileShader(rot3_handle);
#ifdef DEBUG
    printf("---> rot3 Shader:\n");
    debug(rot3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadlfnoise()
{
    int lfnoise_size = strlen(lfnoise_source);
    lfnoise_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(lfnoise_handle, 1, (GLchar **)&lfnoise_source, &lfnoise_size);
    glCompileShader(lfnoise_handle);
#ifdef DEBUG
    printf("---> lfnoise Shader:\n");
    debug(lfnoise_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddbox()
{
    int dbox_size = strlen(dbox_source);
    dbox_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dbox_handle, 1, (GLchar **)&dbox_source, &dbox_size);
    glCompileShader(dbox_handle);
#ifdef DEBUG
    printf("---> dbox Shader:\n");
    debug(dbox_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddbox3()
{
    int dbox3_size = strlen(dbox3_source);
    dbox3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dbox3_handle, 1, (GLchar **)&dbox3_source, &dbox3_size);
    glCompileShader(dbox3_handle);
#ifdef DEBUG
    printf("---> dbox3 Shader:\n");
    debug(dbox3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddvoronoi()
{
    int dvoronoi_size = strlen(dvoronoi_source);
    dvoronoi_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dvoronoi_handle, 1, (GLchar **)&dvoronoi_source, &dvoronoi_size);
    glCompileShader(dvoronoi_handle);
#ifdef DEBUG
    printf("---> dvoronoi Shader:\n");
    debug(dvoronoi_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddquadvoronoi()
{
    int dquadvoronoi_size = strlen(dquadvoronoi_source);
    dquadvoronoi_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dquadvoronoi_handle, 1, (GLchar **)&dquadvoronoi_source, &dquadvoronoi_size);
    glCompileShader(dquadvoronoi_handle);
#ifdef DEBUG
    printf("---> dquadvoronoi Shader:\n");
    debug(dquadvoronoi_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadanalytical_box()
{
    int analytical_box_size = strlen(analytical_box_source);
    analytical_box_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(analytical_box_handle, 1, (GLchar **)&analytical_box_source, &analytical_box_size);
    glCompileShader(analytical_box_handle);
#ifdef DEBUG
    printf("---> analytical_box Shader:\n");
    debug(analytical_box_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadmfnoise()
{
    int mfnoise_size = strlen(mfnoise_source);
    mfnoise_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(mfnoise_handle, 1, (GLchar **)&mfnoise_source, &mfnoise_size);
    glCompileShader(mfnoise_handle);
#ifdef DEBUG
    printf("---> mfnoise Shader:\n");
    debug(mfnoise_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddpolygon()
{
    int dpolygon_size = strlen(dpolygon_source);
    dpolygon_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dpolygon_handle, 1, (GLchar **)&dpolygon_source, &dpolygon_size);
    glCompileShader(dpolygon_handle);
#ifdef DEBUG
    printf("---> dpolygon Shader:\n");
    debug(dpolygon_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddstar()
{
    int dstar_size = strlen(dstar_source);
    dstar_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dstar_handle, 1, (GLchar **)&dstar_source, &dstar_size);
    glCompileShader(dstar_handle);
#ifdef DEBUG
    printf("---> dstar Shader:\n");
    debug(dstar_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddcirclesegment()
{
    int dcirclesegment_size = strlen(dcirclesegment_source);
    dcirclesegment_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dcirclesegment_handle, 1, (GLchar **)&dcirclesegment_source, &dcirclesegment_size);
    glCompileShader(dcirclesegment_handle);
#ifdef DEBUG
    printf("---> dcirclesegment Shader:\n");
    debug(dcirclesegment_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddcircle()
{
    int dcircle_size = strlen(dcircle_source);
    dcircle_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dcircle_handle, 1, (GLchar **)&dcircle_source, &dcircle_size);
    glCompileShader(dcircle_handle);
#ifdef DEBUG
    printf("---> dcircle Shader:\n");
    debug(dcircle_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddlinesegment()
{
    int dlinesegment_size = strlen(dlinesegment_source);
    dlinesegment_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dlinesegment_handle, 1, (GLchar **)&dlinesegment_source, &dlinesegment_size);
    glCompileShader(dlinesegment_handle);
#ifdef DEBUG
    printf("---> dlinesegment Shader:\n");
    debug(dlinesegment_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddlogo210()
{
    int dlogo210_size = strlen(dlogo210_source);
    dlogo210_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dlogo210_handle, 1, (GLchar **)&dlogo210_source, &dlogo210_size);
    glCompileShader(dlogo210_handle);
#ifdef DEBUG
    printf("---> dlogo210 Shader:\n");
    debug(dlogo210_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrand3()
{
    int rand3_size = strlen(rand3_source);
    rand3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rand3_handle, 1, (GLchar **)&rand3_source, &rand3_size);
    glCompileShader(rand3_handle);
#ifdef DEBUG
    printf("---> rand3 Shader:\n");
    debug(rand3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddspline3()
{
    int dspline3_size = strlen(dspline3_source);
    dspline3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dspline3_handle, 1, (GLchar **)&dspline3_source, &dspline3_size);
    glCompileShader(dspline3_handle);
#ifdef DEBUG
    printf("---> dspline3 Shader:\n");
    debug(dspline3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}

void LoadSymbols()
{
    Loadrand();
    updateBar();
    Loadzextrude();
    updateBar();
    Loadstroke();
    updateBar();
    Loadsmoothmin();
    updateBar();
    Loaddhexagonpattern();
    updateBar();
    Loadnormal();
    updateBar();
    Loadrot3();
    updateBar();
    Loadlfnoise();
    updateBar();
    Loaddbox();
    updateBar();
    Loaddbox3();
    updateBar();
    Loaddvoronoi();
    updateBar();
    Loaddquadvoronoi();
    updateBar();
    Loadanalytical_box();
    updateBar();
    Loadmfnoise();
    updateBar();
    Loaddpolygon();
    updateBar();
    Loaddstar();
    updateBar();
    Loaddcirclesegment();
    updateBar();
    Loaddcircle();
    updateBar();
    Loaddlinesegment();
    updateBar();
    Loaddlogo210();
    updateBar();
    Loadrand3();
    updateBar();
    Loaddspline3();
    updateBar();
}
int hexagontunnel_program, hexagontunnel_handle, voronoinet_program, voronoinet_handle, startunnel_program, startunnel_handle, team210_logo_program, team210_logo_handle, broccoli_program, broccoli_handle, boxplane_program, boxplane_handle, doublependulum_program, doublependulum_handle, midikeyboard_program, midikeyboard_handle, wursttunnel_program, wursttunnel_handle;
int hexagontunnel_iTime_location;
hexagontunnel_iFFTWidth_location;
hexagontunnel_iScale_location;
hexagontunnel_iHighScale_location;
hexagontunnel_iNBeats_location;
hexagontunnel_iDial0_location;
hexagontunnel_iDial6_location;
hexagontunnel_iDial7_location;
hexagontunnel_iResolution_location;
hexagontunnel_iFFT_location;
hexagontunnel_iNote_location;
hexagontunnel_iPressure_location;
int voronoinet_iTime_location;
voronoinet_iFFTWidth_location;
voronoinet_iScale_location;
voronoinet_iHighScale_location;
voronoinet_iNBeats_location;
voronoinet_iDial0_location;
voronoinet_iDial6_location;
voronoinet_iDial7_location;
voronoinet_iResolution_location;
voronoinet_iFFT_location;
voronoinet_iNote_location;
voronoinet_iPressure_location;
int startunnel_iTime_location;
startunnel_iFFTWidth_location;
startunnel_iScale_location;
startunnel_iHighScale_location;
startunnel_iNBeats_location;
startunnel_iDial0_location;
startunnel_iDial6_location;
startunnel_iDial7_location;
startunnel_iResolution_location;
startunnel_iFFT_location;
startunnel_iNote_location;
startunnel_iPressure_location;
int team210_logo_iTime_location;
team210_logo_iFFTWidth_location;
team210_logo_iScale_location;
team210_logo_iHighScale_location;
team210_logo_iNBeats_location;
team210_logo_iDial0_location;
team210_logo_iDial6_location;
team210_logo_iDial7_location;
team210_logo_iResolution_location;
team210_logo_iFFT_location;
team210_logo_iNote_location;
team210_logo_iPressure_location;
int broccoli_iTime_location;
broccoli_iFFTWidth_location;
broccoli_iScale_location;
broccoli_iHighScale_location;
broccoli_iNBeats_location;
broccoli_iDial0_location;
broccoli_iDial6_location;
broccoli_iDial7_location;
broccoli_iResolution_location;
broccoli_iFFT_location;
broccoli_iNote_location;
broccoli_iPressure_location;
int boxplane_iTime_location;
boxplane_iFFTWidth_location;
boxplane_iScale_location;
boxplane_iHighScale_location;
boxplane_iNBeats_location;
boxplane_iDial0_location;
boxplane_iDial6_location;
boxplane_iDial7_location;
boxplane_iResolution_location;
boxplane_iFFT_location;
boxplane_iNote_location;
boxplane_iPressure_location;
int doublependulum_iTime_location;
doublependulum_iFFTWidth_location;
doublependulum_iScale_location;
doublependulum_iHighScale_location;
doublependulum_iNBeats_location;
doublependulum_iDial0_location;
doublependulum_iDial6_location;
doublependulum_iDial7_location;
doublependulum_iResolution_location;
doublependulum_iFFT_location;
doublependulum_iNote_location;
doublependulum_iPressure_location;
int midikeyboard_iTime_location;
midikeyboard_iFFTWidth_location;
midikeyboard_iScale_location;
midikeyboard_iHighScale_location;
midikeyboard_iNBeats_location;
midikeyboard_iDial0_location;
midikeyboard_iDial6_location;
midikeyboard_iDial7_location;
midikeyboard_iResolution_location;
midikeyboard_iFFT_location;
midikeyboard_iNote_location;
midikeyboard_iPressure_location;
int wursttunnel_iTime_location;
wursttunnel_iFFTWidth_location;
wursttunnel_iScale_location;
wursttunnel_iHighScale_location;
wursttunnel_iNBeats_location;
wursttunnel_iDial0_location;
wursttunnel_iDial6_location;
wursttunnel_iDial7_location;
wursttunnel_iResolution_location;
wursttunnel_iFFT_location;
wursttunnel_iNote_location;
wursttunnel_iPressure_location;
const int nprograms = 9;

void Loadhexagontunnel()
{
    int hexagontunnel_size = strlen(hexagontunnel_source);
    hexagontunnel_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(hexagontunnel_handle, 1, (GLchar **)&hexagontunnel_source, &hexagontunnel_size);
    glCompileShader(hexagontunnel_handle);
#ifdef DEBUG
    printf("---> hexagontunnel Shader:\n");
    debug(hexagontunnel_handle);
    printf(">>>>\n");
#endif
    hexagontunnel_program = glCreateProgram();
    glAttachShader(hexagontunnel_program,hexagontunnel_handle);
    glAttachShader(hexagontunnel_program,rand_handle);
    glAttachShader(hexagontunnel_program,zextrude_handle);
    glAttachShader(hexagontunnel_program,stroke_handle);
    glAttachShader(hexagontunnel_program,smoothmin_handle);
    glAttachShader(hexagontunnel_program,dhexagonpattern_handle);
    glAttachShader(hexagontunnel_program,normal_handle);
    glAttachShader(hexagontunnel_program,rot3_handle);
    glAttachShader(hexagontunnel_program,lfnoise_handle);
    glLinkProgram(hexagontunnel_program);
#ifdef DEBUG
    printf("---> hexagontunnel Program:\n");
    debugp(hexagontunnel_program);
    printf(">>>>\n");
#endif
    glUseProgram(hexagontunnel_program);
    hexagontunnel_iTime_location = glGetUniformLocation(hexagontunnel_program, "iTime");
    hexagontunnel_iFFTWidth_location = glGetUniformLocation(hexagontunnel_program, "iFFTWidth");
    hexagontunnel_iScale_location = glGetUniformLocation(hexagontunnel_program, "iScale");
    hexagontunnel_iHighScale_location = glGetUniformLocation(hexagontunnel_program, "iHighScale");
    hexagontunnel_iNBeats_location = glGetUniformLocation(hexagontunnel_program, "iNBeats");
    hexagontunnel_iDial0_location = glGetUniformLocation(hexagontunnel_program, "iDial0");
    hexagontunnel_iDial6_location = glGetUniformLocation(hexagontunnel_program, "iDial6");
    hexagontunnel_iDial7_location = glGetUniformLocation(hexagontunnel_program, "iDial7");
    hexagontunnel_iResolution_location = glGetUniformLocation(hexagontunnel_program, "iResolution");
    hexagontunnel_iFFT_location = glGetUniformLocation(hexagontunnel_program, "iFFT");
    hexagontunnel_iNote_location = glGetUniformLocation(hexagontunnel_program, "iNote");
    hexagontunnel_iPressure_location = glGetUniformLocation(hexagontunnel_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loadvoronoinet()
{
    int voronoinet_size = strlen(voronoinet_source);
    voronoinet_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(voronoinet_handle, 1, (GLchar **)&voronoinet_source, &voronoinet_size);
    glCompileShader(voronoinet_handle);
#ifdef DEBUG
    printf("---> voronoinet Shader:\n");
    debug(voronoinet_handle);
    printf(">>>>\n");
#endif
    voronoinet_program = glCreateProgram();
    glAttachShader(voronoinet_program,voronoinet_handle);
    glAttachShader(voronoinet_program,rand_handle);
    glAttachShader(voronoinet_program,zextrude_handle);
    glAttachShader(voronoinet_program,stroke_handle);
    glAttachShader(voronoinet_program,dbox_handle);
    glAttachShader(voronoinet_program,dbox3_handle);
    glAttachShader(voronoinet_program,smoothmin_handle);
    glAttachShader(voronoinet_program,dvoronoi_handle);
    glAttachShader(voronoinet_program,dquadvoronoi_handle);
    glAttachShader(voronoinet_program,analytical_box_handle);
    glAttachShader(voronoinet_program,rot3_handle);
    glAttachShader(voronoinet_program,normal_handle);
    glLinkProgram(voronoinet_program);
#ifdef DEBUG
    printf("---> voronoinet Program:\n");
    debugp(voronoinet_program);
    printf(">>>>\n");
#endif
    glUseProgram(voronoinet_program);
    voronoinet_iTime_location = glGetUniformLocation(voronoinet_program, "iTime");
    voronoinet_iFFTWidth_location = glGetUniformLocation(voronoinet_program, "iFFTWidth");
    voronoinet_iScale_location = glGetUniformLocation(voronoinet_program, "iScale");
    voronoinet_iHighScale_location = glGetUniformLocation(voronoinet_program, "iHighScale");
    voronoinet_iNBeats_location = glGetUniformLocation(voronoinet_program, "iNBeats");
    voronoinet_iDial0_location = glGetUniformLocation(voronoinet_program, "iDial0");
    voronoinet_iDial6_location = glGetUniformLocation(voronoinet_program, "iDial6");
    voronoinet_iDial7_location = glGetUniformLocation(voronoinet_program, "iDial7");
    voronoinet_iResolution_location = glGetUniformLocation(voronoinet_program, "iResolution");
    voronoinet_iFFT_location = glGetUniformLocation(voronoinet_program, "iFFT");
    voronoinet_iNote_location = glGetUniformLocation(voronoinet_program, "iNote");
    voronoinet_iPressure_location = glGetUniformLocation(voronoinet_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loadstartunnel()
{
    int startunnel_size = strlen(startunnel_source);
    startunnel_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(startunnel_handle, 1, (GLchar **)&startunnel_source, &startunnel_size);
    glCompileShader(startunnel_handle);
#ifdef DEBUG
    printf("---> startunnel Shader:\n");
    debug(startunnel_handle);
    printf(">>>>\n");
#endif
    startunnel_program = glCreateProgram();
    glAttachShader(startunnel_program,startunnel_handle);
    glAttachShader(startunnel_program,rand_handle);
    glAttachShader(startunnel_program,lfnoise_handle);
    glAttachShader(startunnel_program,mfnoise_handle);
    glAttachShader(startunnel_program,dbox_handle);
    glAttachShader(startunnel_program,dpolygon_handle);
    glAttachShader(startunnel_program,dstar_handle);
    glAttachShader(startunnel_program,rot3_handle);
    glAttachShader(startunnel_program,stroke_handle);
    glLinkProgram(startunnel_program);
#ifdef DEBUG
    printf("---> startunnel Program:\n");
    debugp(startunnel_program);
    printf(">>>>\n");
#endif
    glUseProgram(startunnel_program);
    startunnel_iTime_location = glGetUniformLocation(startunnel_program, "iTime");
    startunnel_iFFTWidth_location = glGetUniformLocation(startunnel_program, "iFFTWidth");
    startunnel_iScale_location = glGetUniformLocation(startunnel_program, "iScale");
    startunnel_iHighScale_location = glGetUniformLocation(startunnel_program, "iHighScale");
    startunnel_iNBeats_location = glGetUniformLocation(startunnel_program, "iNBeats");
    startunnel_iDial0_location = glGetUniformLocation(startunnel_program, "iDial0");
    startunnel_iDial6_location = glGetUniformLocation(startunnel_program, "iDial6");
    startunnel_iDial7_location = glGetUniformLocation(startunnel_program, "iDial7");
    startunnel_iResolution_location = glGetUniformLocation(startunnel_program, "iResolution");
    startunnel_iFFT_location = glGetUniformLocation(startunnel_program, "iFFT");
    startunnel_iNote_location = glGetUniformLocation(startunnel_program, "iNote");
    startunnel_iPressure_location = glGetUniformLocation(startunnel_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loadteam210_logo()
{
    int team210_logo_size = strlen(team210_logo_source);
    team210_logo_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(team210_logo_handle, 1, (GLchar **)&team210_logo_source, &team210_logo_size);
    glCompileShader(team210_logo_handle);
#ifdef DEBUG
    printf("---> team210_logo Shader:\n");
    debug(team210_logo_handle);
    printf(">>>>\n");
#endif
    team210_logo_program = glCreateProgram();
    glAttachShader(team210_logo_program,team210_logo_handle);
    glAttachShader(team210_logo_program,rand_handle);
    glAttachShader(team210_logo_program,lfnoise_handle);
    glAttachShader(team210_logo_program,mfnoise_handle);
    glAttachShader(team210_logo_program,stroke_handle);
    glAttachShader(team210_logo_program,dcirclesegment_handle);
    glAttachShader(team210_logo_program,dcircle_handle);
    glAttachShader(team210_logo_program,dlinesegment_handle);
    glAttachShader(team210_logo_program,dlogo210_handle);
    glAttachShader(team210_logo_program,zextrude_handle);
    glAttachShader(team210_logo_program,dbox3_handle);
    glAttachShader(team210_logo_program,dhexagonpattern_handle);
    glAttachShader(team210_logo_program,rot3_handle);
    glLinkProgram(team210_logo_program);
#ifdef DEBUG
    printf("---> team210_logo Program:\n");
    debugp(team210_logo_program);
    printf(">>>>\n");
#endif
    glUseProgram(team210_logo_program);
    team210_logo_iTime_location = glGetUniformLocation(team210_logo_program, "iTime");
    team210_logo_iFFTWidth_location = glGetUniformLocation(team210_logo_program, "iFFTWidth");
    team210_logo_iScale_location = glGetUniformLocation(team210_logo_program, "iScale");
    team210_logo_iHighScale_location = glGetUniformLocation(team210_logo_program, "iHighScale");
    team210_logo_iNBeats_location = glGetUniformLocation(team210_logo_program, "iNBeats");
    team210_logo_iDial0_location = glGetUniformLocation(team210_logo_program, "iDial0");
    team210_logo_iDial6_location = glGetUniformLocation(team210_logo_program, "iDial6");
    team210_logo_iDial7_location = glGetUniformLocation(team210_logo_program, "iDial7");
    team210_logo_iResolution_location = glGetUniformLocation(team210_logo_program, "iResolution");
    team210_logo_iFFT_location = glGetUniformLocation(team210_logo_program, "iFFT");
    team210_logo_iNote_location = glGetUniformLocation(team210_logo_program, "iNote");
    team210_logo_iPressure_location = glGetUniformLocation(team210_logo_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loadbroccoli()
{
    int broccoli_size = strlen(broccoli_source);
    broccoli_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(broccoli_handle, 1, (GLchar **)&broccoli_source, &broccoli_size);
    glCompileShader(broccoli_handle);
#ifdef DEBUG
    printf("---> broccoli Shader:\n");
    debug(broccoli_handle);
    printf(">>>>\n");
#endif
    broccoli_program = glCreateProgram();
    glAttachShader(broccoli_program,broccoli_handle);
    glAttachShader(broccoli_program,rand_handle);
    glAttachShader(broccoli_program,lfnoise_handle);
    glAttachShader(broccoli_program,rand3_handle);
    glAttachShader(broccoli_program,zextrude_handle);
    glAttachShader(broccoli_program,stroke_handle);
    glAttachShader(broccoli_program,smoothmin_handle);
    glAttachShader(broccoli_program,dbox3_handle);
    glAttachShader(broccoli_program,normal_handle);
    glLinkProgram(broccoli_program);
#ifdef DEBUG
    printf("---> broccoli Program:\n");
    debugp(broccoli_program);
    printf(">>>>\n");
#endif
    glUseProgram(broccoli_program);
    broccoli_iTime_location = glGetUniformLocation(broccoli_program, "iTime");
    broccoli_iFFTWidth_location = glGetUniformLocation(broccoli_program, "iFFTWidth");
    broccoli_iScale_location = glGetUniformLocation(broccoli_program, "iScale");
    broccoli_iHighScale_location = glGetUniformLocation(broccoli_program, "iHighScale");
    broccoli_iNBeats_location = glGetUniformLocation(broccoli_program, "iNBeats");
    broccoli_iDial0_location = glGetUniformLocation(broccoli_program, "iDial0");
    broccoli_iDial6_location = glGetUniformLocation(broccoli_program, "iDial6");
    broccoli_iDial7_location = glGetUniformLocation(broccoli_program, "iDial7");
    broccoli_iResolution_location = glGetUniformLocation(broccoli_program, "iResolution");
    broccoli_iFFT_location = glGetUniformLocation(broccoli_program, "iFFT");
    broccoli_iNote_location = glGetUniformLocation(broccoli_program, "iNote");
    broccoli_iPressure_location = glGetUniformLocation(broccoli_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loadboxplane()
{
    int boxplane_size = strlen(boxplane_source);
    boxplane_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(boxplane_handle, 1, (GLchar **)&boxplane_source, &boxplane_size);
    glCompileShader(boxplane_handle);
#ifdef DEBUG
    printf("---> boxplane Shader:\n");
    debug(boxplane_handle);
    printf(">>>>\n");
#endif
    boxplane_program = glCreateProgram();
    glAttachShader(boxplane_program,boxplane_handle);
    glLinkProgram(boxplane_program);
#ifdef DEBUG
    printf("---> boxplane Program:\n");
    debugp(boxplane_program);
    printf(">>>>\n");
#endif
    glUseProgram(boxplane_program);
    boxplane_iTime_location = glGetUniformLocation(boxplane_program, "iTime");
    boxplane_iFFTWidth_location = glGetUniformLocation(boxplane_program, "iFFTWidth");
    boxplane_iScale_location = glGetUniformLocation(boxplane_program, "iScale");
    boxplane_iHighScale_location = glGetUniformLocation(boxplane_program, "iHighScale");
    boxplane_iNBeats_location = glGetUniformLocation(boxplane_program, "iNBeats");
    boxplane_iDial0_location = glGetUniformLocation(boxplane_program, "iDial0");
    boxplane_iDial6_location = glGetUniformLocation(boxplane_program, "iDial6");
    boxplane_iDial7_location = glGetUniformLocation(boxplane_program, "iDial7");
    boxplane_iResolution_location = glGetUniformLocation(boxplane_program, "iResolution");
    boxplane_iFFT_location = glGetUniformLocation(boxplane_program, "iFFT");
    boxplane_iNote_location = glGetUniformLocation(boxplane_program, "iNote");
    boxplane_iPressure_location = glGetUniformLocation(boxplane_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loaddoublependulum()
{
    int doublependulum_size = strlen(doublependulum_source);
    doublependulum_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(doublependulum_handle, 1, (GLchar **)&doublependulum_source, &doublependulum_size);
    glCompileShader(doublependulum_handle);
#ifdef DEBUG
    printf("---> doublependulum Shader:\n");
    debug(doublependulum_handle);
    printf(">>>>\n");
#endif
    doublependulum_program = glCreateProgram();
    glAttachShader(doublependulum_program,doublependulum_handle);
    glAttachShader(doublependulum_program,rand_handle);
    glAttachShader(doublependulum_program,lfnoise_handle);
    glAttachShader(doublependulum_program,stroke_handle);
    glAttachShader(doublependulum_program,rot3_handle);
    glAttachShader(doublependulum_program,normal_handle);
    glLinkProgram(doublependulum_program);
#ifdef DEBUG
    printf("---> doublependulum Program:\n");
    debugp(doublependulum_program);
    printf(">>>>\n");
#endif
    glUseProgram(doublependulum_program);
    doublependulum_iTime_location = glGetUniformLocation(doublependulum_program, "iTime");
    doublependulum_iFFTWidth_location = glGetUniformLocation(doublependulum_program, "iFFTWidth");
    doublependulum_iScale_location = glGetUniformLocation(doublependulum_program, "iScale");
    doublependulum_iHighScale_location = glGetUniformLocation(doublependulum_program, "iHighScale");
    doublependulum_iNBeats_location = glGetUniformLocation(doublependulum_program, "iNBeats");
    doublependulum_iDial0_location = glGetUniformLocation(doublependulum_program, "iDial0");
    doublependulum_iDial6_location = glGetUniformLocation(doublependulum_program, "iDial6");
    doublependulum_iDial7_location = glGetUniformLocation(doublependulum_program, "iDial7");
    doublependulum_iResolution_location = glGetUniformLocation(doublependulum_program, "iResolution");
    doublependulum_iFFT_location = glGetUniformLocation(doublependulum_program, "iFFT");
    doublependulum_iNote_location = glGetUniformLocation(doublependulum_program, "iNote");
    doublependulum_iPressure_location = glGetUniformLocation(doublependulum_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loadmidikeyboard()
{
    int midikeyboard_size = strlen(midikeyboard_source);
    midikeyboard_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(midikeyboard_handle, 1, (GLchar **)&midikeyboard_source, &midikeyboard_size);
    glCompileShader(midikeyboard_handle);
#ifdef DEBUG
    printf("---> midikeyboard Shader:\n");
    debug(midikeyboard_handle);
    printf(">>>>\n");
#endif
    midikeyboard_program = glCreateProgram();
    glAttachShader(midikeyboard_program,midikeyboard_handle);
    glLinkProgram(midikeyboard_program);
#ifdef DEBUG
    printf("---> midikeyboard Program:\n");
    debugp(midikeyboard_program);
    printf(">>>>\n");
#endif
    glUseProgram(midikeyboard_program);
    midikeyboard_iTime_location = glGetUniformLocation(midikeyboard_program, "iTime");
    midikeyboard_iFFTWidth_location = glGetUniformLocation(midikeyboard_program, "iFFTWidth");
    midikeyboard_iScale_location = glGetUniformLocation(midikeyboard_program, "iScale");
    midikeyboard_iHighScale_location = glGetUniformLocation(midikeyboard_program, "iHighScale");
    midikeyboard_iNBeats_location = glGetUniformLocation(midikeyboard_program, "iNBeats");
    midikeyboard_iDial0_location = glGetUniformLocation(midikeyboard_program, "iDial0");
    midikeyboard_iDial6_location = glGetUniformLocation(midikeyboard_program, "iDial6");
    midikeyboard_iDial7_location = glGetUniformLocation(midikeyboard_program, "iDial7");
    midikeyboard_iResolution_location = glGetUniformLocation(midikeyboard_program, "iResolution");
    midikeyboard_iFFT_location = glGetUniformLocation(midikeyboard_program, "iFFT");
    midikeyboard_iNote_location = glGetUniformLocation(midikeyboard_program, "iNote");
    midikeyboard_iPressure_location = glGetUniformLocation(midikeyboard_program, "iPressure");
    progress += .2/(float)nprograms;
}

void Loadwursttunnel()
{
    int wursttunnel_size = strlen(wursttunnel_source);
    wursttunnel_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(wursttunnel_handle, 1, (GLchar **)&wursttunnel_source, &wursttunnel_size);
    glCompileShader(wursttunnel_handle);
#ifdef DEBUG
    printf("---> wursttunnel Shader:\n");
    debug(wursttunnel_handle);
    printf(">>>>\n");
#endif
    wursttunnel_program = glCreateProgram();
    glAttachShader(wursttunnel_program,wursttunnel_handle);
    glAttachShader(wursttunnel_program,rand_handle);
    glAttachShader(wursttunnel_program,lfnoise_handle);
    glAttachShader(wursttunnel_program,stroke_handle);
    glAttachShader(wursttunnel_program,dspline3_handle);
    glAttachShader(wursttunnel_program,rot3_handle);
    glAttachShader(wursttunnel_program,normal_handle);
    glLinkProgram(wursttunnel_program);
#ifdef DEBUG
    printf("---> wursttunnel Program:\n");
    debugp(wursttunnel_program);
    printf(">>>>\n");
#endif
    glUseProgram(wursttunnel_program);
    wursttunnel_iTime_location = glGetUniformLocation(wursttunnel_program, "iTime");
    wursttunnel_iFFTWidth_location = glGetUniformLocation(wursttunnel_program, "iFFTWidth");
    wursttunnel_iScale_location = glGetUniformLocation(wursttunnel_program, "iScale");
    wursttunnel_iHighScale_location = glGetUniformLocation(wursttunnel_program, "iHighScale");
    wursttunnel_iNBeats_location = glGetUniformLocation(wursttunnel_program, "iNBeats");
    wursttunnel_iDial0_location = glGetUniformLocation(wursttunnel_program, "iDial0");
    wursttunnel_iDial6_location = glGetUniformLocation(wursttunnel_program, "iDial6");
    wursttunnel_iDial7_location = glGetUniformLocation(wursttunnel_program, "iDial7");
    wursttunnel_iResolution_location = glGetUniformLocation(wursttunnel_program, "iResolution");
    wursttunnel_iFFT_location = glGetUniformLocation(wursttunnel_program, "iFFT");
    wursttunnel_iNote_location = glGetUniformLocation(wursttunnel_program, "iNote");
    wursttunnel_iPressure_location = glGetUniformLocation(wursttunnel_program, "iPressure");
    progress += .2/(float)nprograms;
}

void LoadPrograms()
{
    Loadhexagontunnel();
    updateBar();
    Loadvoronoinet();
    updateBar();
    Loadstartunnel();
    updateBar();
    Loadteam210_logo();
    updateBar();
    Loadbroccoli();
    updateBar();
    Loadboxplane();
    updateBar();
    Loaddoublependulum();
    updateBar();
    Loadmidikeyboard();
    updateBar();
    Loadwursttunnel();
    updateBar();
}
#endif
