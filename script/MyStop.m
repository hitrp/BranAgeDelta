print('-dpdf',OutFile);
system(sprintf('pdfcrop %s.pdf %s.pdf',OutFile,OutFile));
system(sprintf('open %s.pdf',OutFile));
