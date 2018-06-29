/**通过外部程序打开Excel文件**/


public class StartExcel{
	public static void start( )
	{
		
		String cmd = "cmd /C D:\\POI_test.xls";
		Process p;
		try {
			p = Runtime.getRuntime().exec(cmd);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
