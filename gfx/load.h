/* Generated with shader-compressor by NR4/Team210. */
#ifndef LOAD_H
#define LOAD_H
const char * load_frag =
"/* Tunguska by Team210 - 64k Intro at Solskogen 2k19\n"
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
"#version 130\n"
"\n"
"uniform float iTime, iProgress;\n"
"uniform vec2 iResolution;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0, ry = 1.0;\n"
"\n"
"// Hash function\n"
"void rand(in vec2 x, out float num)\n"
"{\n"
"    x += 400.;\n"
"    num = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\n"
"// 2D box\n"
"void dbox2(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"	vec2 da = abs(x)-b;\n"
"	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\n"
"// Simple numbers\n"
"void dnumber(in vec2 x, in int number, in vec2 size, out float d)\n"
"{\n"
"    if(number == 0)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"        dbox2(x, vec2(.4,.6)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 1)\n"
"    {\n"
"        dbox2(x,vec2(.4,1.)*size, d);\n"
"    }\n"
"    else if(number == 2)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"        dbox2(x+vec2(.3,-.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    	dbox2(x-vec2(.3,-.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 3)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"        dbox2(x+vec2(.3,-.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    	dbox2(x+vec2(.3,.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 4)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"        dbox2(x-vec2(0.,.6)*size, vec2(.5,.5)*size, da);\n"
"        d = max(d, -da);\n"
"    	dbox2(x-vec2(-.3,-.7)*size, vec2(.8,.4)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 5)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"        dbox2(x+vec2(.3,.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    	dbox2(x-vec2(.3,.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 6)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"        dbox2(x+vec2(0.,.4)*size, vec2(.5,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    	dbox2(x-vec2(.3,.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 7)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"    	dbox2(x-vec2(-.3,-.3)*size, vec2(.8,.9)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 8)\n"
"    {\n"
"        dbox2(x, size, d);\n"
"        float da;\n"
"        dbox2(x+vec2(0.,.4)*size, vec2(.5,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    	dbox2(x-vec2(0.,.4)*size, vec2(.5,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"    else if(number == 9)\n"
"    {\n"
"        dbox2(-x, size, d);\n"
"        float da;\n"
"        dbox2(-x+vec2(0.,.4)*size, vec2(.5,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    	dbox2(-x-vec2(.3,.4)*size, vec2(.8,.2)*size, da);\n"
"        d = max(d, -da);\n"
"    }\n"
"}\n"
"\n"
"// Convert progress to %.2d\n"
"void dprogress(in vec2 x, in vec2 size, out float d)\n"
"{\n"
"    int n = int(floor(10.*clamp(iProgress,0.,.9999)));\n"
"    dnumber(x, n, size, d);\n"
"    float da;\n"
"    n = int(floor(100.*clamp(iProgress,0.,.9999))) - 10 * n;\n"
"    dnumber(x-2.1*size.x*c.xy, n, size, da);\n"
"    d = min(d,da);\n"
"}\n"
"\n"
"void dteam210(in vec2 x, in float size, out float d)\n"
"{\n"
"    dbox2(x, vec2(.2,1.)*size, d);\n"
"    d = min(d, abs(length(x-vec2(.8,0.)*size)-.8*size)-.2*size);\n"
"    d = min(d, abs(length(x+vec2(.8,0.)*size)-.8*size)-.2*size);\n"
"    float da;\n"
"    dbox2(x+2.85*size*c.xy,2.*vec2(size), da);\n"
"    d = max(d,-da);\n"
"}\n"
"\n"
"// compute distance to regular triangle\n"
"void dtriangle2(in vec2 uv, in float r, out float d)\n"
"{\n"
"    float dp = 2.*pi/3.;\n"
"    vec2 p0 = r*vec2(cos(pi/2.), -sin(pi/2.)),\n"
"        p1 = r*vec2(cos(pi/2.+dp), -sin(pi/2.+dp)),\n"
"        p2 = r*vec2(cos(pi/2.+2.*dp), -sin(pi/2.+2.*dp)), \n"
"        pd = p2-p1;\n"
"    \n"
"    d = min(dot(uv-p0,c.xz*(p1-p0).yx),dot(uv-p1, pd.yx*c.xz));\n"
"	d = min(d, dot(uv-p2, (p0-p2).yx*c.xz))/length(pd);\n"
"}\n"
"\n"
"// distance to gear\n"
"void dgear(in vec2 x, in vec2 r, in float n, out float d)\n"
"{\n"
"    float p = atan(x.y,x.x);\n"
"    p = mod(p, 2.*pi/n)*n/2./pi;\n"
"    d = mix(length(x)-r.x, length(x)-r.y, step(p,.5));\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    ry = 1.5/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    // Progress text\n"
"    float d, v;\n"
"    dprogress(uv+.51*c.xy, vec2(.03,.04), d);\n"
"    col = mix(col,vec3(0.75,0.20,0.26), smoothstep(ry, -ry, d));\n"
"    \n"
"    // Bar outline\n"
"    dbox2(uv, vec2(.4,.04), d);\n"
"    d = abs(d)-.002;\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry, -ry, d));\n"
"    \n"
"    // Bar content\n"
"    dbox2(uv+(.42-.42*clamp(iProgress,0.,.9999))*c.xy, vec2(clamp(iProgress,0.,.9999)*.42,.06), d);\n"
"    vec3 fc = vec3(0.76,0.20,0.25);\n"
"    col = mix(col, fc, smoothstep(ry, -ry, d+.03));\n"
"    \n"
"    // 210 Logo\n"
"    dteam210(uv-.13*c.yx+.025*c.xy-(-.42+.84*clamp(iProgress,0.,.9999))*c.xy, .05, d);\n"
"    col = mix(col, vec3(0.24,0.24,0.24), smoothstep(ry,-ry,d));\n"
"    d = abs(d-.01)-.002;\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry,d));\n"
"    \n"
"    // Red progress triangle\n"
"    dtriangle2(uv-.055*c.yx+.025*c.xy-(-.42+.84*clamp(iProgress,0.,.9999))*c.xy, .01, d);\n"
"    col = mix(col, vec3(0.76,0.20,0.25), smoothstep(ry,-ry, -d));\n"
"    \n"
"    // Upper phone task bar background\n"
"    dbox2(uv-.475*c.yx, vec2(a,.025), d);\n"
"    col = mix(col, vec3(0.24,0.24,0.24), smoothstep(ry,-ry, d));\n"
"    \n"
"    // Attention sign\n"
"    dtriangle2((uv-vec2(-.48*a,.47))*c.xz, .025, d);\n"
"    d = -d;\n"
"    dbox2(uv-vec2(-.48*a,.477), .5*vec2(.005,.017), v);\n"
"    d = max(d,-v);\n"
"    dbox2(uv-vec2(-.48*a,.462), .5*vec2(.005, .005), v);\n"
"    d = max(d,-v);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry,d));\n"
"    \n"
"    // Battery Block\n"
"    dbox2((uv-vec2(.45*a,.475))*c.xz, vec2(.035, .018), d);\n"
"    dbox2((uv-vec2(.47*a,.475))*c.xz, vec2(.008, .01), v);\n"
"    d = min(d,v);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry,d));\n"
"    dbox2((uv-vec2(.439*a,.475))*c.xz, .9*vec2(.015, .018), d);\n"
"    col = mix(col, vec3(0.76,0.20,0.25), smoothstep(ry,-ry,d));\n"
"    \n"
"    // Network information\n"
"    dbox2((uv-vec2(.39*a,.466))*c.xz, vec2(.0055, .009), d);\n"
"    dbox2((uv-vec2(.398*a,.469))*c.xz, vec2(.0055, .012), v);\n"
"	d = min(d,v);\n"
"    dbox2((uv-vec2(.406*a,.472))*c.xz, vec2(.0055, .015), v);\n"
"	d = min(d,v);\n"
"    dbox2((uv-vec2(.414*a,.475))*c.xz, vec2(.0055, .018), v);\n"
"	d = min(d,v);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry,d));\n"
"    \n"
"    // Wifi\n"
"    d = length(uv-vec2(.36*a,.46))-.0025;\n"
"    d = min(d, abs(length(uv-vec2(.36*a,.46))-.01)-.003);\n"
"    d = min(d, abs(length(uv-vec2(.36*a,.46))-.02)-.003);\n"
"    d = min(d, abs(length(uv-vec2(.36*a,.46))-.03)-.003);\n"
"    mat2 m = mat2(cos(pi/4.),-sin(pi/4.), sin(pi/4.), cos(pi/4.));\n"
"    dbox2(m*(uv-vec2(.36*a,.5)), 2.*vec2(.015), v);\n"
"    d = max(d,v);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry,d));\n"
"    \n"
"    // Playback bar\n"
"    dbox2(uv-.41*c.yx*c.xz, vec2(.47*a,.003), d);\n"
"    col = mix(col, c.xyy, smoothstep(ry,-ry, d));\n"
"    dbox2(uv-.41*c.yx*c.xz-.4*c.xy, vec2(.47*a-.4,.003), d);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry, d));\n"
"    \n"
"    // Play symbol\n"
"    dtriangle2(((uv-vec2(-.44*a,-.455))).yx*c.xz, .025, d);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry, -d));\n"
"    \n"
"    // Next symbol\n"
"    dtriangle2(((uv-vec2(-.4*a,-.455))).yx*c.xz, .025, d);\n"
"    dbox2(uv-vec2(-.385*a,-.455), vec2(.003, .024), v);\n"
"    d = min(-d,v);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry, d));\n"
"    \n"
"    // Mute symbol\n"
"    dtriangle2(((uv-vec2(-.35*a,-.455))).yx, .025, d);\n"
"    dbox2(uv-vec2(-.36*a,-.455), vec2(.007, .007), v);\n"
"    d = min(-d,v);\n"
"    d = min(d, length(uv-vec2(-.345*a,-.455))-.011);\n"
"    dbox2(uv-vec2(-.344*a,-.455), vec2(.003, .024), v);\n"
"    d = max(d,-v);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry, d));\n"
"    \n"
"    // Settings\n"
"    dgear(uv-vec2(.4*a,-.455), vec2(.02,.025), 8., d);\n"
"    d = max(d, -length(uv-vec2(.4*a,-.455))+.01);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry, d));\n"
"    \n"
"    // Full Screen\n"
"    vec2 y = mod(uv, .03)-.015;\n"
"    dbox2(y, vec2(.01), d);\n"
"    d = abs(d)-.003;\n"
"    dbox2(uv-vec2(.44*a-.0025,-.45), vec2(.02), v);\n"
"	d = max(d,v);\n"
"    col = mix(col, vec3(0.92,0.89,0.84), smoothstep(ry,-ry, d));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.0,1.0),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\n"
;
#endif
