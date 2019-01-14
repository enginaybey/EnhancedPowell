function[xmin]=goldsec(fname,a,b,delta)
t=0.5*(sqrt(5.0)-1);
c=b-t*(b-a);
d=a+t*(b-a);
fc=feval(fname,c);
fd=feval(fname,d);
while abs(b-a)>delta
  if fc>fd
    a=c;
    c=d;
    fc=fd;
    d=a+t*(b-a);
    fd=feval(fname,d);
  else
    b=d;
    d=c;
    fd=fc;
    c=b-t*(b-a);
    fc=feval(fname,c);
  end
end
xmin=(c+d)/2;

