//Generated with Symbolize (c) 2019 Alexander Kraus <nr4@z10.info>.
#ifndef SYMBOLIZE_H
#define SYMBOLIZE_H

int rand_handle, zextrude_handle, stroke_handle, smoothmin_handle, dhexagonpattern_handle, normal_handle, rot3_handle, lfnoise_handle, dbox_handle, dbox3_handle, dvoronoi_handle, dquadvoronoi_handle, analytical_box_handle;
const int nsymbols = 13;
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
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
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
"    x.z -= 1.*iTime;\n"
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
"            rot3(na*1.e3*vec3(1.1,1.5,1.9),RR);\n"
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
"uniform vec2 iResolution;\n"
"uniform sampler1D iFFT;\n"
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
"    d.x = max(d.x, dbound);\n"
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
"    float d = 0.;// -(o.z-.12)/dir.z;\n"
"    int N = 250,\n"
"        i;\n"
"    \n"
"    // Graph\n"
"    analytical_box(o,dir,vec3(.3*c.xx,.2),d);\n"
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
"            \n"
"            \n"
"            mat3 RR;\n"
"            rot3(na*1.e3*vec3(1.1,1.5,1.9),RR);\n"
"            col = mix(mix(.0,.3,clamp(x.z/.3,0.,1.))*(.5+.5*mat)*c.xxx,(1.+.8*mat)*abs(RR*RR*vec3(.7,.5,.26)),step(x.z,.08));\n"
"            col = mix(col,(1.+.8*mat)*abs(RR*vec3(.6,.12,.06)),step(.19,x.z));\n"
"\n"
"            col = mix((.5+.5*mat)*col,(1.+.8*mat)*abs(RR*vec3(0.89,0.44,0.23)),(.5+.5*sin(x.z))*step(.19,x.z));\n"
"            col = mix(col,vec3(0.25,0.23,0.21),(.5+.5*cos(4.*x.z+mat))*step(.19,x.z));\n"
"            \n"
"            //dvoronoi(mix(1.,4.,na)*1.01*vec2(pi,3.)*vec2(phi,x.z-iTime),dhex,ind);\n"
"            //stroke(dhex, .3, dhex);\n"
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
}
int hexagontunnel_program, hexagontunnel_handle, voronoinet_program, voronoinet_handle;
int hexagontunnel_iTime_location;
hexagontunnel_iFFTWidth_location;
hexagontunnel_iScale_location;
hexagontunnel_iHighScale_location;
hexagontunnel_iNBeats_location;
hexagontunnel_iResolution_location;
hexagontunnel_iFFT_location;
int voronoinet_iTime_location;
voronoinet_iFFTWidth_location;
voronoinet_iScale_location;
voronoinet_iHighScale_location;
voronoinet_iNBeats_location;
voronoinet_iResolution_location;
voronoinet_iFFT_location;
const int nprograms = 2;

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
    hexagontunnel_iResolution_location = glGetUniformLocation(hexagontunnel_program, "iResolution");
    hexagontunnel_iFFT_location = glGetUniformLocation(hexagontunnel_program, "iFFT");
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
    voronoinet_iResolution_location = glGetUniformLocation(voronoinet_program, "iResolution");
    voronoinet_iFFT_location = glGetUniformLocation(voronoinet_program, "iFFT");
    progress += .2/(float)nprograms;
}

void LoadPrograms()
{
    Loadhexagontunnel();
    updateBar();
    Loadvoronoinet();
    updateBar();
}
#endif
