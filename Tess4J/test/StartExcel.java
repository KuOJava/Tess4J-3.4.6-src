import java.io.File;

/**ͨ���ⲿ�����Excel�ļ�**/
import java.io.IOException;

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
