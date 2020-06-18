model A
    class Pin
		Real v;
		flow Real i;	
    end Pin;
    Pin a[5, 5, 5],b[5, 5, 5],c[5],r;
	class BB
		Real tt;
	end BB;
	BB qw;
equation
	for i in 1:3 loop
      for j in 1:5 loop
		connect(a[i, j],b[i, j]);	
      end for;
	end for;

	for i in 1:5 loop
		connect(b[i],c[i ]);	
	end for;

    connect(a[1:3, 1, 3], b[1:5, 1, 3]);

	for i in 1:5 loop
		connect(c[i],r);	
	end for;	
	a[1].v = 10;
end A;

