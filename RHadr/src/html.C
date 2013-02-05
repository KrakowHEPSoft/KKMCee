{
    gSystem.Load("libKorw.so");
    THtml html;
    html.SetSourceDir("./");
    html.SetOutputDir("./dok/");
    html.MakeAll();
}
