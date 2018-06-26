import java.io.File;
/***创建Excel文件，并将数据存入文件中***/
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.io.FileUtils;
import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;

import jxl.write.Label;

public class PoiExcel {
public void poi(ArrayList<String>temp1,ArrayList<String>temp2)
{
	File file=new File("D:\\POI_test.xls");
	String[] title= {"企业名称","企业注册号"};
	HSSFWorkbook workbook=new HSSFWorkbook ();
	HSSFSheet sheet=workbook.createSheet();
	HSSFRow row=sheet.createRow(0); 
	HSSFCell cell=null;
	//表头
	for(int i=0;i<title.length;i++)
	{
		cell=row.createCell(i);
		cell.setCellValue(title[i]);
	}
	
	for(int i=1;i<=temp1.size();i++)
	{

		HSSFRow nextrow=sheet.createRow(i);
		HSSFCell cell2=nextrow.createCell(0);
		cell2.setCellValue(temp1.get(i-1));
		
		cell2=nextrow.createCell(1);
		cell2.setCellValue(temp2.get(i-1));
	}
	

	try {
		file.createNewFile();
		FileOutputStream stream=FileUtils.openOutputStream(file);
		workbook.write(stream);
		stream.close();
	
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}

}
}