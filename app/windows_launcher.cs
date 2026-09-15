// Small launcher: uses the installed Python and the selected current worktree.
using System;
using System.Diagnostics;
using System.IO;
using System.Windows.Forms;
class AlFatigueLauncher {
    static string Quote(string value) { return ((char)34)+value+((char)34); }
    [STAThread] static int Main(string[] args) {
        try {
            string[] paths=File.ReadAllLines(Path.Combine(AppDomain.CurrentDomain.BaseDirectory,"launcher.paths"));
            if(paths.Length!=2 || !File.Exists(paths[0]) || !File.Exists(Path.Combine(paths[1],"app","desktop_ui.py"))) throw new Exception("UI installation paths are unavailable.");
            bool check=args.Length==1 && args[0]=="--check";
            if(args.Length>1 || (args.Length==1 && !check && !args[0].EndsWith(".ftgsim",StringComparison.OrdinalIgnoreCase))) throw new Exception("Open an .ftgsim project.");
            string command="-m app.desktop_ui"+(check?" --smoke":args.Length==1?" "+Quote(Path.GetFullPath(args[0])):"");
            ProcessStartInfo start=new ProcessStartInfo(paths[0],command);
            start.WorkingDirectory=paths[1]; start.UseShellExecute=false; start.CreateNoWindow=true;
            using(Process child=Process.Start(start)) { if(check) { child.WaitForExit(); return child.ExitCode; } }
            return 0;
        } catch(Exception e) { MessageBox.Show(e.Message,"Al Fatigue UI"); return 1; }
    }
}
