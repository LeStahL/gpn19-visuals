//Generated with Symbolize (c) 2019 Alexander Kraus <nr4@z10.info>.
#ifndef SYMBOLIZE_H
#define SYMBOLIZE_H

int rand_handle, zextrude_handle, stroke_handle, smoothmin_handle, dhexagonpattern_handle, normal_handle, rot3_handle, lfnoise_handle, mfnoise_handle, dvoronoi_handle, add_handle, dbox_handle, line_handle;
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
const char *add_source = "void add(in vec2 sda, in vec2 sdb, out vec2 sdf)\n"
"{\n"
"    sdf = mix(sda, sdb, step(sdb.x, sda.x));\n"
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
const char *line_source = "// Distance to line segment\n"
"void line(in vec3 x, in vec3 p1, in vec3 p2, out float dst)\n"
"{\n"
"    vec3 d = p2-p1;\n"
"    dst = length(x-mix(p1, p2, clamp(dot(x-p1, d)/dot(d,d),0.,1.)));\n"
"}\n"
"\0";
const char *decayingfactory_source = "/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
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
"\n"
"float mat;\n"
"void scene(in vec3 x, out vec2 d)\n"
"{\n"
"    d = c.xx;\n"
"    \n"
"    x.z -= .5*iTime;\n"
"    \n"
"    float phi = atan(x.y, x.x),\n"
"        dhex;\n"
"    vec2 ind;\n"
"    dhexagonpattern(2.*1.01*vec2(pi,3.)*vec2(phi,x.z),dhex,ind);\n"
"    rand(ind,mat);\n"
"    stroke(dhex, .1, dhex);\n"
"    mat *= (1.-iScale);\n"
"    //d.x = length(x.xy)-mix(.5,.45,smoothstep(.1,-.1,dhex));\n"
"    d.x = min(d.x, length(x.xy)+.2*mat-mix(.5,.55+.2*mat,smoothstep(.2,-.2,dhex)));\n"
"    \n"
"    stroke(d.x, .04, d.x);\n"
"}\n"
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
"            col = mix((.5+.5*mat)*c.xxx,(1.+.8*mat)*vec3(0.89,0.44,0.23),.5+.5*sin(x.z));\n"
"            col = mix(col,vec3(0.25,0.23,0.21),.5+.5*cos(4.*x.z+mat));\n"
"//             mat3 RR;\n"
"//             rot3(.4*vec3(1.1,1.5,1.9)*iTime-.1*x.z, RR);\n"
"//             col = abs(RR*col);\n"
"        }\n"
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
const char *fogforest_source = "/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
" * Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
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
"// Hash function\n"
"void rand(in vec2 x, out float num);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void dbox(in vec2 p, in vec2 b, out float dst);\n"
"void line(in vec3 x, in vec3 p1, in vec3 p2, out float dst);\n"
"void stroke(in float d, in float s, out float dst);\n"
"\n"
"// Hash function\n"
"void rand3(in vec3 x, out float num)\n"
"{\n"
"    num = fract(sin(dot(sign(x)*abs(x) ,vec3(12.9898,78.233,45.1232)))*43758.5453);\n"
"}\n"
"\n"
"// Arbitrary-frequency 2D noise\n"
"void lfnoise3(in vec3 t, out float num)\n"
"{\n"
"    vec3 i = floor(t);\n"
"    t = fract(t);\n"
"    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise\n"
"    t = smoothstep(c.yyy, c.xxx, t); // TODO: add this for faster value noise\n"
"    vec2 v1, v2, v3, v4;\n"
"    rand3(i, v1.x);\n"
"    rand3(i+c.xyy, v1.y);\n"
"    rand3(i+c.yxy, v2.x);\n"
"    rand3(i+c.xxy, v2.y);\n"
"    rand3(i+c.yyx, v3.x);\n"
"    rand3(i+c.xyx, v3.y);\n"
"    rand3(i+c.yxx, v4.x);\n"
"    rand3(i+c.xxx, v4.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    v3 = c.zz+2.*mix(v3, v4, t.y);\n"
"    v2.x = -1.+2.*mix(v1.x, v1.y, t.x);\n"
"    v2.y = -1.+2.*mix(v3.x, v3.y, t.x);\n"
"    num = mix(v2.x, v2.y, t.z);\n"
"}\n"
"\n"
"// Make noise multi-frequency\n"
"void mfnoise3(in vec3 x, in float fmin, in float fmax, in float alpha, out float dst)\n"
"{\n"
"    dst = 0.;\n"
"    float a = 1., nf = 0.;\n"
"    for(float f = fmin; f<fmax; f = f*2.)\n"
"    {\n"
"        float buf;\n"
"        lfnoise3(f*x, buf);\n"
"        dst += a*buf;\n"
"        a *= alpha;\n"
"        nf += 1.;\n"
"    }\n"
"    dst *= (1.-alpha)/(1.-pow(alpha, nf));\n"
"}\n"
"\n"
"// Scene\n"
"void scene(in vec3 x, out vec2 d)\n"
"{\n"
"    x.y += .1*iTime;\n"
"    float n;\n"
"    \n"
"    // Floor\n"
"    mfnoise(x.xy, 4.,4.e2,.35, n);\n"
"    d = vec2(x.z-.05*n,1.);\n"
"    \n"
"    // Trees\n"
"    float v, da;\n"
"    vec2 vi;\n"
"    dvoronoi(x.xy, v, vi);\n"
"    vec2 r;\n"
"    rand(vi, r.x);\n"
"    rand(vi+1301., r.y);\n"
"    vec2 y = x.xy-vi, \n"
"        n2;\n"
"    lfnoise(3.*x.z*c.xx-r,n2.x);\n"
"    lfnoise(4.*x.z*c.xx+33.1*r, n2.y);\n"
"    da = length(y-.01*n2)-.07*mix(1.,.7,smoothstep(0., 1.,clamp(x.z*3.,0.,1.)));\n"
"    //zextrude(x.z,-da,10.,da);\n"
"    //add(d, vec2(da,2.), d);\n"
"    smoothmin(d.x,da,.2,d.x);\n"
"    d.y = mix(1.,2.,step(.1*n,x.z));\n"
"    \n"
"    // smaller branches\n"
"    float z = mod(x.z,.05)-.025, zi = (x.z-z)/.05;\n"
"    vec2 rp;//= vec2(0;\n"
"    rand(zi*c.xx+r,rp.x);\n"
"    rand(zi*c.xx+r+1332.,rp.y);\n"
"    rp *= vec2(1.,2.*pi);\n"
"\n"
"    float nz;\n"
"    lfnoise(5.*length(y-.01*n2)*c.xx-33.*zi, nz);\n"
"    \n"
"\n"
"    line(vec3(y-.01*n2, z+.01*nz), c.yyy, vec3(rp.x*vec2(cos(rp.y),sin(rp.y)),.05*rp.x), da);\n"
"    stroke(da, mix(1.,.1,smoothstep(0.,1.,clamp(length(vec3(y-.01*n2, z))/.7,0.,1.)))*.01*(.3+rp.x+n2.x), da);\n"
"    smoothmin(d.x,da,.05,d.x);\n"
"}\n"
"\n"
"// Normal\n"
"const float dx = 5.e-4;\n"
"void normal(in vec3 x, out vec3 n);\n"
"\n"
"// Texture\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"  x.y += .1*iTime;\n"
"  float n;\n"
"    mfnoise(x.xy, 4.e0,4.e2,.65, n);\n"
"  col = mix(vec3(0.55,0.69,0.37), vec3(0.00,0.02,0.04), .9+.1*n);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"     // Set up coordinates\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
" \n"
"    if(length(uv) > .5)\n"
"    {\n"
"        fragColor = vec4(col, 0.);\n"
"        return;\n"
"    }\n"
" \n"
"    // Camera setup\n"
"    float pp = .3*iTime;\n"
"    vec3 o = c.yzy+.2*c.yyx, \n"
"        t = c.yyy+.3*c.yyx,\n"
"        dir = normalize(t-o),\n"
"        r = normalize(c.xyy),\n"
"        u = normalize(cross(r,dir)),\n"
"        n,\n"
"        x;\n"
"    t += uv.x*r + uv.y*u;\n"
"    dir = normalize(t-o);\n"
"    vec2 s;\n"
"    float d = 0.;//-(o.z-.05)/dir.z;\n"
"    int N = 350,\n"
"        i;\n"
"    \n"
"    // Raymarching\n"
"    for(i=0; i<N; ++i)\n"
"    {\n"
"        x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4*max(d*d,1.)) break;\n"
"        if(d>10.)break;\n"
"        d += min(.01*max(d,1.),s.x);\n"
"    }\n"
"    \n"
"    // Illumination\n"
"    vec3 l = normalize(x+c.yxx);\n"
"    if(i<N)\n"
"    {\n"
"	    normal(x,n);\n"
"        //colorize(x.xy, col);\n"
"    }\n"
"\n"
"    \n"
"    if(s.y == 2.)//Treess\n"
"    {\n"
"      \n"
"    col = .2*vec3(0.05,0.12,0.12)\n"
"        + .2*vec3(0.05,0.12,0.12)*abs(dot(l,n))\n"
"        + .6*vec3(0.04,0.13,0.12)*abs(pow(dot(reflect(-l,n),dir),3.));\n"
"    }\n"
"    if(s.y == 1.)\n"
"    {\n"
"        colorize(x.xy,col);\n"
"            .5*col\n"
"            + .2*col*abs(dot(l,n))\n"
"            +.6*col*abs(pow(dot(reflect(-l,n),dir),3.));\n"
"    }\n"
"    vec3 c1 =  mix(vec3(0.91,0.87,0.68),vec3(0.07,0.21,0.21),clamp(length(uv),0.,1.));\n"
"    float noiz;\n"
"    mfnoise3(x,1.,100.,.65,noiz);\n"
"    noiz = .5+.5*noiz;\n"
"    //noiz *= smoothstep(.3,.5,clamp(x.z,0.,1.));\n"
"    c1 = mix(c1, vec3(0.29,0.60,0.47), noiz);\n"
"    col = mix(col, c1, clamp(d/10.,0.,1.));\n"
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
void Loadadd()
{
    int add_size = strlen(add_source);
    add_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(add_handle, 1, (GLchar **)&add_source, &add_size);
    glCompileShader(add_handle);
#ifdef DEBUG
    printf("---> add Shader:\n");
    debug(add_handle);
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
void Loadline()
{
    int line_size = strlen(line_source);
    line_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(line_handle, 1, (GLchar **)&line_source, &line_size);
    glCompileShader(line_handle);
#ifdef DEBUG
    printf("---> line Shader:\n");
    debug(line_handle);
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
    Loadmfnoise();
    updateBar();
    Loaddvoronoi();
    updateBar();
    Loadadd();
    updateBar();
    Loaddbox();
    updateBar();
    Loadline();
    updateBar();
}
int decayingfactory_program, decayingfactory_handle, fogforest_program, fogforest_handle;
int decayingfactory_iTime_location;
decayingfactory_iFFTWidth_location;
decayingfactory_iScale_location;
decayingfactory_iHighScale_location;
decayingfactory_iNBeats_location;
decayingfactory_iResolution_location;
decayingfactory_iFFT_location;
int fogforest_iTime_location;
fogforest_iFFTWidth_location;
fogforest_iScale_location;
fogforest_iHighScale_location;
fogforest_iNBeats_location;
fogforest_iResolution_location;
fogforest_iFFT_location;
const int nprograms = 2;

void Loaddecayingfactory()
{
    int decayingfactory_size = strlen(decayingfactory_source);
    decayingfactory_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(decayingfactory_handle, 1, (GLchar **)&decayingfactory_source, &decayingfactory_size);
    glCompileShader(decayingfactory_handle);
#ifdef DEBUG
    printf("---> decayingfactory Shader:\n");
    debug(decayingfactory_handle);
    printf(">>>>\n");
#endif
    decayingfactory_program = glCreateProgram();
    glAttachShader(decayingfactory_program,decayingfactory_handle);
    glAttachShader(decayingfactory_program,rand_handle);
    glAttachShader(decayingfactory_program,zextrude_handle);
    glAttachShader(decayingfactory_program,stroke_handle);
    glAttachShader(decayingfactory_program,smoothmin_handle);
    glAttachShader(decayingfactory_program,dhexagonpattern_handle);
    glAttachShader(decayingfactory_program,normal_handle);
    glAttachShader(decayingfactory_program,rot3_handle);
    glLinkProgram(decayingfactory_program);
#ifdef DEBUG
    printf("---> decayingfactory Program:\n");
    debugp(decayingfactory_program);
    printf(">>>>\n");
#endif
    glUseProgram(decayingfactory_program);
    decayingfactory_iTime_location = glGetUniformLocation(decayingfactory_program, "iTime");
    decayingfactory_iFFTWidth_location = glGetUniformLocation(decayingfactory_program, "iFFTWidth");
    decayingfactory_iScale_location = glGetUniformLocation(decayingfactory_program, "iScale");
    decayingfactory_iHighScale_location = glGetUniformLocation(decayingfactory_program, "iHighScale");
    decayingfactory_iNBeats_location = glGetUniformLocation(decayingfactory_program, "iNBeats");
    decayingfactory_iResolution_location = glGetUniformLocation(decayingfactory_program, "iResolution");
    decayingfactory_iFFT_location = glGetUniformLocation(decayingfactory_program, "iFFT");
    progress += .2/(float)nprograms;
}

void Loadfogforest()
{
    int fogforest_size = strlen(fogforest_source);
    fogforest_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fogforest_handle, 1, (GLchar **)&fogforest_source, &fogforest_size);
    glCompileShader(fogforest_handle);
#ifdef DEBUG
    printf("---> fogforest Shader:\n");
    debug(fogforest_handle);
    printf(">>>>\n");
#endif
    fogforest_program = glCreateProgram();
    glAttachShader(fogforest_program,fogforest_handle);
    glAttachShader(fogforest_program,rand_handle);
    glAttachShader(fogforest_program,lfnoise_handle);
    glAttachShader(fogforest_program,mfnoise_handle);
    glAttachShader(fogforest_program,dvoronoi_handle);
    glAttachShader(fogforest_program,smoothmin_handle);
    glAttachShader(fogforest_program,add_handle);
    glAttachShader(fogforest_program,zextrude_handle);
    glAttachShader(fogforest_program,dbox_handle);
    glAttachShader(fogforest_program,line_handle);
    glAttachShader(fogforest_program,stroke_handle);
    glAttachShader(fogforest_program,normal_handle);
    glLinkProgram(fogforest_program);
#ifdef DEBUG
    printf("---> fogforest Program:\n");
    debugp(fogforest_program);
    printf(">>>>\n");
#endif
    glUseProgram(fogforest_program);
    fogforest_iTime_location = glGetUniformLocation(fogforest_program, "iTime");
    fogforest_iFFTWidth_location = glGetUniformLocation(fogforest_program, "iFFTWidth");
    fogforest_iScale_location = glGetUniformLocation(fogforest_program, "iScale");
    fogforest_iHighScale_location = glGetUniformLocation(fogforest_program, "iHighScale");
    fogforest_iNBeats_location = glGetUniformLocation(fogforest_program, "iNBeats");
    fogforest_iResolution_location = glGetUniformLocation(fogforest_program, "iResolution");
    fogforest_iFFT_location = glGetUniformLocation(fogforest_program, "iFFT");
    progress += .2/(float)nprograms;
}

void LoadPrograms()
{
    Loaddecayingfactory();
    updateBar();
    Loadfogforest();
    updateBar();
}
#endif
