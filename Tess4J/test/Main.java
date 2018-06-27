import java.awt.Component;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

//������������ͼ���û������д������
public class Main {
	JTextField inputPath;
	private JButton choosePath;
	private JButton confirm;
	private JLabel situation;
	public static JTextArea showMessage;
	private JLabel jumpToExcel;
	private JButton clickToCopy;

	public Main() {
		// ����һ�������Ϊ��������
		// JFrame.setDefaultLookAndFeelDecorated(true);
		JFrame f = new JFrame("����ͼƬʶ��ϵͳ");
		// f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		// ��ȡ����е��������
		Container contentPane = f.getContentPane();
		// �������񲼾�
		contentPane.setLayout(new GridBagLayout());
		// ʹ�����񲼾�Լ����
		GridBagConstraints c = new GridBagConstraints();
		// �������,���������ĵ����ṩ��ͼ�ν����
		// '����·��'�ı��򣺴��ڵ�1�У�ռ3��;'ѡ��·��'��ť�����ڵ�1�У���4�У�'ȷ��'��ť�����ڵ�1�У���5��
		// 'ʶ�����'��ʾ�򣺴��ڵ�2�У�ռ������
		// ��ʾ��Ϣ���ı��򣺴��ڵ�3,4�У�ռ������
		// '������������Excel�ĸ�Ŀ¼'�����ӣ����ڵ�5�У�ռ3�У�'��������ļ�'��ť�����ڵ�5�У�ռ����

		// '����·��'���ı���
		inputPath = new JTextField(1);
		// ��һ�е�һ�п�ʼ
		c.gridx = 0;
		c.gridy = 0;
		// ��3��
		c.gridwidth = 3;
		c.gridheight = 1;
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 3;
		c.weighty = 1;
		// ���õ���¼�

		// ����c��������������
		contentPane.add(inputPath, c);

		// 'ѡ��·��'��ť
		choosePath = new JButton("ѡ��·��");
		// ѡ��·������¼�
		choosePath.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				Component parent = null;
				int returnVal = chooser.showOpenDialog(parent);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					inputPath.setText(chooser.getSelectedFile().getAbsolutePath());
				}
			}
		});
		// ���ڵ�һ�е�����
		c.gridx = 3;
		c.gridy = 0;
		c.gridwidth = 1;
		c.weightx = 1;
		c.weighty = 1;
		// c.gridheight=1;----���ø�ֵ
		contentPane.add(choosePath, c);
		// 'ȷ��'��ť
		confirm = new JButton("ȷ��");
		// ȷ����ť�ĵ���¼�
		confirm.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// ������ȡ·����ȡ����
					String path = inputPath.getText();
					if (path.equals("")) {
						showMessage.setText("��û�������κ�·��\n");
					} else {
						showMessage.setText("����ʶ�𡣡���������");
						new Thread(new Runnable() {

							@Override
							public void run() {
								// TODO Auto-generated method stub
								String message="";
								try {
									message = Tess.find(path, Tess.depth);
									
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
								PoiExcel poiexcel = new PoiExcel();
								poiexcel.poi(Tess.temp1, Tess.temp2);// ����Excel���
								showMessage.setText(message);
							}
						}).start();

					}
			}
		});
		c.gridx = 4;
		// c.gridy=4;----���ø�ֵ
		// c.gridwidth=1;----���ø�ֵ
		// c.gridheight=1;----���ø�ֵ
		c.weightx = 1;
		c.weighty = 1;
		contentPane.add(confirm, c);

		// 'ʶ�����:'��ʾ��
		situation = new JLabel("ʶ�����:");
		c.gridx = 0;
		c.gridy = 1;
		c.weightx = 1;
		c.weighty = 1;
		// c.gridwidth=1;----���ø�ֵ
		// c.gridheight=1;----���ø�ֵ
		contentPane.add(situation, c);

		// '��ʾ��Ϣ�Ŀ�'
		showMessage = new JTextArea("");
		showMessage.setSize(200, 200);
		// c.gridx=0;
		c.gridy = 2;
		c.gridwidth = 5;
		c.gridheight = 2;
		c.weightx = 5;
		c.weighty = 2;
		contentPane.add(showMessage, c);

		// '������������Excel�ĸ�Ŀ¼'
		jumpToExcel = new JLabel("������������Excel�ĸ�Ŀ¼");
		c.gridx = 0;
		c.gridy = 4;
		c.gridwidth = 3;
		c.gridheight = 1;
		c.weightx = 3;
		c.weighty = 1;
		contentPane.add(jumpToExcel, c);

		// '��������ļ�'
		clickToCopy = new JButton("��������ļ�");
		c.gridx = 3;
		c.gridwidth = 2;
		c.weightx = 2;
		c.weighty = 1;
		contentPane.add(clickToCopy, c);
		// ��ʾ����
		f.pack();
		f.setSize(800, 1000);
		f.setVisible(true);
	}

	public static void main(String[] args) {
		new Main();
	}
}
