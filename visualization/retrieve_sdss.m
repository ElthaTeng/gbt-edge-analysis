pk=parseCsvTable('galaxy_list.csv');
name=pk.galname;
disp(length(name));

for i=1:length(name)
    disp(name(i));
    v=sprintf('https://nsatlas.org/getAtlas.html?search=name&name=%s&radius=10.0&submit_form=Submit',name{i});
    data=webread(v);
    ix=strfind(data,'<img src=');
    if ~isempty(ix)
        imgv=sscanf(data(ix(1):end),'<img src="%s" />');
        img=webread(imgv(1:(end-1)));
        imwrite(img,['cutouts/' name{i} '.jpg'],'JPEG');
        clf;
        imagesc(img); axis image; drawnow;
    end
end